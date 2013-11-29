# ======================================================================
# Atomistica - Interatomic potential library
# https://github.com/pastewka/atomistica
# Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
# See the AUTHORS file in the top-level Atomistica directory.
#
# Copyright (2005-2013) Fraunhofer IWM
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
            mods += [ l.split(':')[0:4] ]
        l = f.readline()
    f.close()

    return mods

###

def write_factory_f90(mods, str, fn):
    f = open(fn, "w")

    f.write("#include \"macros.inc\"\n\n" +
            "module %s_factory\n" % str +
            '  use libAtoms_module\n' +
            '  use particles\n' +
            '  use neighbors\n')
    for f90name, f90class, name, register_data_ex in mods:
        f.write('  use %s\n' % f90name)
    f.write("  implicit none\n\n" +
            "contains\n\n")

    for f90name, f90class, name, register_data_ex in mods:
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
                "  if (.not. associated(this_fptr))   stop '[lammps_%s_free] *this_fptr* is NULL.'\n" % f90name +
                "  call del(this_fptr)\n" +
                "  deallocate(this_fptr)\n"  +
                "endsubroutine lammps_%s_free\n\n\n" % f90name)

        if register_data_ex == "1":
            raise RuntimeError('MDCORE potential {0} has register_data '
                               'interface. Cannot interface that to LAMMPS.'
                               .format(f90name))

        f.write("subroutine lammps_%s_init_without_parameters(this_cptr) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[%s_init_without_parameters] *this_fptr* is NULL.'\n" % f90name +
                "  call init(this_fptr)\n" +
                "endsubroutine lammps_%s_init_without_parameters\n\n\n" % f90name)

        f.write("subroutine lammps_%s_del(this_cptr) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[%s_del] *this_fptr* is NULL.'\n" % f90name +
                "  call del(this_fptr)\n" +
                "endsubroutine lammps_%s_del\n\n\n" % f90name)

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

        f.write("subroutine lammps_%s_energy_and_forces(this_cptr, p_cptr, nl_cptr, epot, f, wpot, epot_per_at_cptr, epot_per_bond_cptr, f_per_bond_cptr, wpot_per_at_cptr, wpot_per_bond_cptr, error) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(c_ptr), value :: p_cptr\n" +
                "  type(c_ptr), value :: nl_cptr\n" +
                "  real(c_double), intent(out) :: epot\n" +
                "  real(c_double) :: f(3, *)\n" +
                "  real(c_double) :: wpot(3, 3)\n" +
                "  type(c_ptr), value :: epot_per_at_cptr\n" +
                "  type(c_ptr), value :: epot_per_bond_cptr\n" +
                "  type(c_ptr), value :: f_per_bond_cptr\n" +
                "  type(c_ptr), value :: wpot_per_at_cptr\n" +
                "  type(c_ptr), value :: wpot_per_bond_cptr\n" +
                "  integer(c_int), intent(out) :: error\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  type(particles_t), pointer :: p\n" +
                "  type(neighbors_t), pointer :: nl\n" +
                "  real(c_double), pointer :: epot_per_at(:), epot_per_bond(:), f_per_bond(:,:)\n" +
                "  real(c_double), pointer :: wpot_per_at(:,:), wpot_per_bond(:,:)\n" +
                "  error = ERROR_NONE\n" +
                "  nullify(epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond)\n" +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  call c_f_pointer(p_cptr, p)\n" +
                "  call c_f_pointer(nl_cptr, nl)\n" +
                "  if (c_associated(epot_per_at_cptr)) then\n" +
                "     call c_f_pointer(epot_per_at_cptr, epot_per_at, [p%nat])\n" +
                "  endif\n" +
                "  if (c_associated(epot_per_bond_cptr)) then\n" +
                "     call c_f_pointer(epot_per_bond_cptr, epot_per_bond, [nl%neighbors_size])\n" +
                "  endif\n" +
                "  if (c_associated(f_per_bond_cptr)) then\n" +
                "     call c_f_pointer(f_per_bond_cptr, f_per_bond, [3,p%nat])\n" +
                "  endif\n" +
                "  if (c_associated(wpot_per_at_cptr)) then\n" +
                "     call c_f_pointer(wpot_per_at_cptr, wpot_per_at, [6,p%nat])\n" +
                "  endif\n" +
                "  if (c_associated(wpot_per_bond_cptr)) then\n" +
                "     call c_f_pointer(wpot_per_bond_cptr, wpot_per_bond, [6,nl%neighbors_size])\n" +
                "  endif\n" +
                "  if (.not. associated(this_fptr))   stop '[lammps_%s_energy_and_forces] *this_fptr* is NULL.'\n" % f90name +
                "  if (.not. associated(p))   stop '[lammps_%s_energy_and_forces] *p* is NULL.'\n" % f90name +
                "  if (.not. associated(nl))   stop '[lammps_%s_energy_and_forces] *nl* is NULL.'\n" % f90name +
                "  call energy_and_forces(this_fptr, p, nl, epot, f, wpot, &\n" +
                "    epot_per_at=epot_per_at, epot_per_bond=epot_per_bond, f_per_bond=f_per_bond, &\n" +
                "    wpot_per_at=wpot_per_at, wpot_per_bond=wpot_per_bond, ierror=error)\n" +
                "endsubroutine lammps_%s_energy_and_forces\n\n\n" % f90name)
    
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
    for f90name, f90class, name, register_data_ex in mods:
        s += "void lammps_%s_new(void **, section_t *, section_t **);\n" % f90name
        s += "void lammps_%s_free(void *);\n" % f90name
        s += "void lammps_%s_init_without_parameters(void *);\n" % f90name
        s += "void lammps_%s_del(void *);\n" % f90name
        s += "void lammps_%s_bind_to(void *, void *, void *, int *);\n" % f90name
        s += "void lammps_%s_energy_and_forces(void *, void *, void *, double *, double *, double *, double *, double *, double *, double *, double *, int *);\n" % f90name

    d["prototypes"] = s

    #
    # Classes
    #

    s = "%s_class_t %s_classes[N_POTENTIAL_CLASSES] = {\n" % ( str, str )
    for f90name, f90class, name, register_data_ex in mods:
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
