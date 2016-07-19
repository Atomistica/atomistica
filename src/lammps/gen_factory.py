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

import itertools
import sys

###

def file_from_template(templatefn, fn, keywords):
    ftmp = open(templatefn, "r")
    fout = open(fn, "w")

    l = ftmp.readline()
    while l:
        fout.write(l % keywords)
        l = ftmp.readline()

    ftmp.close()
    fout.close()

###

def read_module_list(fn):
    mods = [ ]

    f = open(fn, "r")
    l = f.readline()
    while l:
        l  = l.strip()
        if len(l) > 0 and l[0] != '!' and l[0] != '#':
            s = l.split(':')[0:5]
            s[4] = s[4].split(',')
            mods += [s]
        l = f.readline()
    f.close()

    return mods

###

def switch_optargs(funcstr, optargs):
    s = ''
    if len(optargs) == 0:
        s += '  call %s\n' % (funcstr % '')
    else:
        for perm in itertools.product(*([[True,False]]*len(optargs))):
            cond = '.true.'
            args = ''
            for condp, arg in zip(perm, optargs):
                if condp:
                    cond += ' .and. associated(%s)' % arg
                    args += '%s=%s, ' % (arg, arg)
                else:
                    cond += ' .and. .not. associated(%s)' % arg
            s += '  if (%s) then\n' % cond
            s += '     call %s\n' % (funcstr % args)
            s += '  else\n'
        s += '     stop "Fatal internal error: Dispatch should not have ended up here."\n'
        for perm in itertools.product(*([[True,False]]*len(optargs))):
            s += '  endif\n'
    return s

###

def write_factory_f90(mods, str, fn):
    f = open(fn, "w")

    f.write("#include \"macros.inc\"\n\n" +
            "module %s_factory\n" % str +
            '  use libAtoms_module\n' +
            '  use particles\n' +
            '  use neighbors\n')
    for f90name, f90class, name, features, methods in mods:
        f.write('  use %s\n' % f90name)
    f.write("  implicit none\n\n" +
            "contains\n\n")

    for f90name, f90class, name, features, methods in mods:
        features = set(features.split(','))
        f.write("subroutine lammps_%s_new(this_cptr, cfg, m) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), intent(out) :: this_cptr\n" +
                "  type(c_ptr), value :: cfg\n" +
                "  type(c_ptr), intent(out) :: m\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  allocate(this_fptr)\n"  +
                "  if (.not. associated(this_fptr))   stop '[lammps_%s_new] *this_fptr* is NULL.'\n" % f90name +
                "  call register(this_fptr, cfg, m)\n" +
                "  this_cptr = c_loc(this_fptr)\n" +
                "endsubroutine lammps_%s_new\n\n\n" % f90name)


        f.write("subroutine lammps_%s_free(this_cptr) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[lammps_%s_free] *this_fptr* is NULL.'\n" % f90name)
        if 'del' in methods:
            f.write("  call del(this_fptr)\n")
        f.write("#ifndef __bg__\n" +
                "  deallocate(this_fptr)\n"  +
                "#endif\n" +
                "endsubroutine lammps_%s_free\n\n\n" % f90name)

        if 'register_data' in methods:
            raise RuntimeError('MDCORE potential {0} has register_data '
                               'interface. Cannot interface that to LAMMPS.'
                               .format(f90name))

        f.write("subroutine lammps_%s_init_without_parameters(this_cptr) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[%s_init_without_parameters] *this_fptr* is NULL.'\n" % f90name)
        if 'init' in methods:
            f.write("  call init(this_fptr)\n")
        f.write("endsubroutine lammps_%s_init_without_parameters\n\n\n" % f90name)

        f.write("subroutine lammps_%s_del(this_cptr) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[%s_del] *this_fptr* is NULL.'\n" % f90name)
        if 'del' in methods:
            f.write("  call del(this_fptr)\n")
        f.write("endsubroutine lammps_%s_del\n\n\n" % f90name)

        f.write("subroutine lammps_%s_bind_to(this_cptr, p_cptr, nl_cptr, error) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(c_ptr), value :: p_cptr\n" +
                "  type(c_ptr), value :: nl_cptr\n" +
                "  integer(c_int), intent(out) :: error\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  type(particles_t), pointer :: p\n" +
                "  type(neighbors_t), pointer :: nl\n" +
                "  error = ERROR_NONE\n" +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  call c_f_pointer(p_cptr, p)\n" +
                "  call c_f_pointer(nl_cptr, nl)\n" +
                "  if (.not. associated(this_fptr))   stop '[lammps_%s_bind_to] *this_fptr* is NULL.'\n" % f90name +
                "  if (.not. associated(p))   stop '[lammps_%s_bind_to] *p* is NULL.'\n" % f90name +
                "  if (.not. associated(nl))   stop '[lammps_%s_bind_to] *nl* is NULL.'\n" % f90name +
                "  call bind_to(this_fptr, p, nl, ierror=error)\n" +
                "endsubroutine lammps_%s_bind_to\n\n\n" % f90name)

        f.write("subroutine lammps_%s_energy_and_forces(this_cptr, p_cptr, nl_cptr, epot, f, wpot, mask_cptr, epot_per_at_cptr, wpot_per_at_cptr, error) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(c_ptr), value :: p_cptr\n" +
                "  type(c_ptr), value :: nl_cptr\n" +
                "  real(c_double), intent(out) :: epot\n" +
                "  real(c_double) :: f(3, *)\n" +
                "  real(c_double) :: wpot(3, 3)\n" +
                "  type(c_ptr), value :: mask_cptr\n" +
                "  type(c_ptr), value :: epot_per_at_cptr\n" +
                "  type(c_ptr), value :: wpot_per_at_cptr\n" +
                "  integer(c_int), intent(out) :: error\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  type(particles_t), pointer :: p\n" +
                "  type(neighbors_t), pointer :: nl\n")
        if 'mask' in features:
                f.write("  integer(c_int), pointer :: mask(:)\n")
        if 'per_at' in features:
                f.write("  real(c_double), pointer :: epot_per_at(:), wpot_per_at(:,:)\n")
        f.write("  error = ERROR_NONE\n" +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  call c_f_pointer(p_cptr, p)\n" +
                "  call c_f_pointer(nl_cptr, nl)\n" +
                "  if (.not. associated(this_fptr))   stop '[lammps_%s_energy_and_forces] *this_fptr* is NULL.'\n" % f90name +
                "  if (.not. associated(p))   stop '[lammps_%s_energy_and_forces] *p* is NULL.'\n" % f90name +
                "  if (.not. associated(nl))   stop '[lammps_%s_energy_and_forces] *nl* is NULL.'\n" % f90name)
        optargs = []
        if 'mask' in features:
            f.write("  call c_f_pointer(mask_cptr, mask, [p%nat])\n")
            optargs += ['mask']
        else:
            f.write('  if (c_associated(mask_cptr)) then\n' +
                    '     RETURN_ERROR("*mask* argument present but not supported by potential %s.", error)\n' % name +
                    '  endif\n')
        if 'per_at' in features:
            f.write("  call c_f_pointer(epot_per_at_cptr, epot_per_at, [p%nat])\n" +
                    "  call c_f_pointer(wpot_per_at_cptr, wpot_per_at, [6,p%nat])\n")
            optargs += ['epot_per_at', 'wpot_per_at']
        else:
            f.write('  if (c_associated(epot_per_at_cptr)) then\n' +
                    '     RETURN_ERROR("*epot_per_at* argument present but not supported by potential %s.", error)\n' % name +
                    '  endif\n' +
                    '  if (c_associated(wpot_per_at_cptr)) then\n' +
                    '     RETURN_ERROR("*wpot_per_at* argument present but not supported by potential %s.", error)\n' % name +
                    '  endif\n')
        f.write(switch_optargs("energy_and_forces(this_fptr, p, nl, epot, f, wpot, %sierror=error)", optargs))
        f.write("endsubroutine lammps_%s_energy_and_forces\n\n\n" % f90name)
    
    f.write("endmodule %s_factory\n" % str)
    f.close()

###

def write_factory_c(mods, str, c_dispatch_template, c_dispatch_file,
                    h_dispatch_template, h_dispatch_file):

    d = { }

    d["disclaimer"] = "This file has been autogenerated. DO NOT MODIFY."
    d["name"] = str
    d["n_classes"] = len(mods)

    #
    # Prototypes
    #

    s = ""
    for f90name, f90class, name, features, methods in mods:
        s += "void lammps_%s_new(void **, section_t *, section_t **);\n" % f90name
        s += "void lammps_%s_free(void *);\n" % f90name
        s += "void lammps_%s_init_without_parameters(void *);\n" % f90name
        s += "void lammps_%s_del(void *);\n" % f90name
        s += "void lammps_%s_bind_to(void *, void *, void *, int *);\n" % f90name
        s += "void lammps_%s_energy_and_forces(void *, void *, void *, double *, double *, double *, int *, double *, double *, int *);\n" % f90name

    d["prototypes"] = s

    #
    # Classes
    #

    s = "%s_class_t %s_classes[N_POTENTIAL_CLASSES] = {\n" % ( str, str )
    for f90name, f90class, name, features, methods in mods:
        s += "  {\n"
        s += "    \"%s\",\n" % name
        s += "    lammps_%s_new,\n" % f90name
        s += "    lammps_%s_free,\n" % f90name
        s += "    lammps_%s_init_without_parameters,\n" % f90name
        s += "    lammps_%s_del,\n" % f90name
        s += "    lammps_%s_bind_to,\n" % f90name
        s += "    lammps_%s_energy_and_forces,\n" % f90name
        s += "  },\n"
            
    s = s[:-2] + "\n};\n"

    d["classes"] = s  

    #
    # Write the dispatch module
    #

    d["dispatch_header"] = h_dispatch_file.split('/')[-1]

    file_from_template(c_dispatch_template, c_dispatch_file, d)
    file_from_template(h_dispatch_template, h_dispatch_file, d)

###

if __name__ == '__main__':
    srcdir, compiler, machine, system  = sys.argv[1:5]

    mods = read_module_list("potentials.classes")
    write_factory_f90(mods, "potential", "potentials_factory_f90.f90")
    write_factory_c(mods, "potential", 
                    srcdir + "/factory.template.c", "potentials_factory_c.c",
                    srcdir + "/factory.template.h", "potentials_factory_c.h")
