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
            mods += [ l.split(':')[0:5] ]
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
        f.write("subroutine python_%s_new(this_cptr, cfg, m) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), intent(out) :: this_cptr\n" +
                "  type(c_ptr), value :: cfg\n" +
                "  type(c_ptr), intent(out) :: m\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  allocate(this_fptr)\n"  +
                "  if (.not. associated(this_fptr))   stop '[python_%s_new] *this_fptr* is NULL.'\n" % f90name +
                "  call register(this_fptr, cfg, m)\n" +
                "  this_cptr = c_loc(this_fptr)\n" +
                "endsubroutine python_%s_new\n\n\n" % f90name)

        f.write("subroutine python_%s_free(this_cptr) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[python_%s_free] *this_fptr* is NULL.'\n" % f90name)
        if 'del' in methods:
            f.write("  call del(this_fptr)\n")
        f.write("  deallocate(this_fptr)\n"  +
                "endsubroutine python_%s_free\n\n\n" % f90name)

        if 'register_data' in methods:
            f.write("subroutine python_%s_register_data(this_cptr, p_cptr, error) bind(C)\n" % f90name +
                    "  use, intrinsic :: iso_c_binding\n\n" +
                    "  implicit none\n\n" +
                    "  type(c_ptr), value :: this_cptr\n" +
                    "  type(c_ptr), value :: p_cptr\n" +
                    "  integer(c_int), intent(out) :: error\n\n" +
                    "  type(%s_t), pointer :: this_fptr\n" % f90name +
                    "  type(particles_t), pointer :: p\n" +
                    "  error = ERROR_NONE\n" +
                    "  call c_f_pointer(this_cptr, this_fptr)\n" +
                    "  call c_f_pointer(p_cptr, p)\n" +
                    "  if (.not. associated(this_fptr))   stop '[python_%s_register_data] *this_fptr* is NULL.'\n" % f90name +
                    "  if (.not. associated(p))   stop '[python_%s_register_data] *p* is NULL.'\n" % f90name +
                    "  call register_data(this_fptr, p, ierror=error)\n" +
                    "endsubroutine python_%s_register_data\n\n\n" % f90name)
        else:
            f.write("subroutine python_%s_register_data(this_cptr, p_cptr, error) bind(C)\n" % f90name +
                    "  use, intrinsic :: iso_c_binding\n\n" +
                    "  implicit none\n\n" +
                    "  type(c_ptr), value :: this_cptr\n" +
                    "  type(c_ptr), value :: p_cptr\n" +
                    "  integer(c_int), intent(out) :: error\n\n" +
                    "  type(%s_t), pointer :: this_fptr\n" % f90name +
                    "  type(particles_t), pointer :: p\n" +
                    "  error = ERROR_NONE\n" +
                    "  call c_f_pointer(this_cptr, this_fptr)\n" +
                    "  call c_f_pointer(p_cptr, p)\n" +
                    "  if (.not. associated(this_fptr))   stop '[python_%s_register_data] *this_fptr* is NULL.'\n" % f90name +
                    "  if (.not. associated(p))   stop '[python_%s_register_data] *p* is NULL.'\n" % f90name +
                    "endsubroutine python_%s_register_data\n\n\n" % f90name)

        f.write("subroutine python_%s_init_without_parameters(this_cptr, error) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  integer(c_int), intent(out) :: error\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  error=ERROR_NONE\n" +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[%s_init_without_parameters] *this_fptr* is NULL.'\n" % f90name)
        if 'init' in methods:
            f.write("  call init(this_fptr)\n")
        f.write("endsubroutine python_%s_init_without_parameters\n\n\n" % f90name)

        if 'set_coulomb' in methods:
            f.write("subroutine python_%s_set_Coulomb(this_cptr, coul_cptr, error) bind(C)\n" % f90name +
                    "  use, intrinsic :: iso_c_binding\n\n" +
                    "  implicit none\n\n" +
                    "  type(c_ptr), value :: this_cptr\n" +
                    "  type(c_ptr), value :: coul_cptr\n" +
                    "  integer(c_int), intent(out) :: error\n\n" +
                    "  type(%s_t), pointer :: this_fptr\n" % f90name +
                    "  error = ERROR_NONE\n" +
                    "  call c_f_pointer(this_cptr, this_fptr)\n" +
                    "  if (.not. associated(this_fptr))   stop '[python_%s_set_Coulomb] *this_fptr* is NULL.'\n" % f90name +
                    "  call set_Coulomb(this_fptr, coul_cptr, ierror=error)\n" +
                    "endsubroutine python_%s_set_Coulomb\n\n\n" % f90name)

        if 'get_dict' in methods:
            f.write("subroutine python_%s_get_dict(this_cptr, dict_cptr, error) bind(C)\n" % f90name +
                    "  use, intrinsic :: iso_c_binding\n\n" +
                    "  implicit none\n\n" +
                    "  type(c_ptr), value :: this_cptr\n" +
                    "  type(c_ptr), value :: dict_cptr\n" +
                    "  integer(c_int), intent(out) :: error\n\n" +
                    "  type(%s_t), pointer :: this_fptr\n" % f90name +
                    "  type(ptrdict_t) :: dict\n" +
                    "  error = ERROR_NONE\n" +
                    "  call c_f_pointer(this_cptr, this_fptr)\n" +
                    "  dict%ptrdict = dict_cptr\n" +
                    "  if (.not. associated(this_fptr))   stop '[python_%s_get_dict] *this_fptr* is NULL.'\n" % f90name +
                    "  call get_dict(this_fptr, dict, error=error)\n" +
                    "endsubroutine python_%s_get_dict\n\n\n" % f90name)

        if 'get_per_bond_property' in methods:
            f.write("subroutine python_%s_get_per_bond_property(this_cptr, p_cptr, nl_cptr, propstr, propout, error) bind(C)\n" % f90name +
                    "  use, intrinsic :: iso_c_binding\n\n" +
                    "  implicit none\n\n" +
                    "  type(c_ptr), value :: this_cptr\n" +
                    "  type(c_ptr), value :: p_cptr\n" +
                    "  type(c_ptr), value :: nl_cptr\n" +
                    "  type(c_ptr), value :: propstr\n" +
                    "  real(c_double), intent(out) :: propout(*)\n\n" +
                    "  integer(c_int), intent(out) :: error\n\n" +
                    "  type(%s_t), pointer :: this_fptr\n" % f90name +
                    "  type(particles_t), pointer :: p\n" +
                    "  type(neighbors_t), pointer :: nl\n" +
                    "  error = ERROR_NONE\n" +
                    "  call c_f_pointer(this_cptr, this_fptr)\n" +
                    "  call c_f_pointer(p_cptr, p)\n" +
                    "  call c_f_pointer(nl_cptr, nl)\n" +
                    "  if (.not. associated(this_fptr))   stop '[python_%s_get_per_bond_property] *this_fptr* is NULL.'\n" % f90name +
                    "  call get_per_bond_property(this_fptr, p, nl, a2s(c_f_string(propstr)), propout, error=error)\n" +
                    "endsubroutine python_%s_get_per_bond_property\n\n\n" % f90name)

        f.write("subroutine python_%s_bind_to(this_cptr, p_cptr, nl_cptr, error) bind(C)\n" % f90name +
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
                "  if (.not. associated(this_fptr))   stop '[python_%s_bind_to] *this_fptr* is NULL.'\n" % f90name +
                "  if (.not. associated(p))   stop '[python_%s_bind_to] *p* is NULL.'\n" % f90name +
                "  if (.not. associated(nl))   stop '[python_%s_bind_to] *nl* is NULL.'\n" % f90name +
                "  call bind_to(this_fptr, p, nl, ierror=error)\n" +
                "endsubroutine python_%s_bind_to\n\n\n" % f90name)

        s = """
subroutine python_%s_energy_and_forces(this_cptr, p_cptr, nl_cptr, &
  q, epot, f, wpot, mask_cptr, epot_per_at_cptr, epot_per_bond_cptr, &
  f_per_bond_cptr, wpot_per_at_cptr, wpot_per_bond_cptr, error) bind(C)
  use, intrinsic :: iso_c_binding

  implicit none

  type(c_ptr), value :: this_cptr
  type(c_ptr), value :: p_cptr
  type(c_ptr), value :: nl_cptr
  real(c_double) :: q(*)
  real(c_double), intent(out) :: epot
  real(c_double) :: f(3, *)
  real(c_double) :: wpot(3, 3)
  type(c_ptr), value :: mask_cptr
  type(c_ptr), value :: epot_per_at_cptr
  type(c_ptr), value :: epot_per_bond_cptr
  type(c_ptr), value :: f_per_bond_cptr
  type(c_ptr), value :: wpot_per_at_cptr
  type(c_ptr), value :: wpot_per_bond_cptr
  integer(c_int), intent(out) :: error

  type(%s_t), pointer :: this_fptr
  type(particles_t), pointer :: p
  type(neighbors_t), pointer :: nl
        """ % ( f90name, f90name )
        if 'mask' in features:
            s += """
  integer(c_int), pointer :: mask(:)
            """
        if 'per_at' in features:
            s += """
  real(c_double), pointer :: epot_per_at(:)
  real(c_double), pointer :: wpot_per_at(:, :, :)
            """
        if 'per_bond' in features:
            s += """
  real(c_double), pointer :: epot_per_bond(:)
  real(c_double), pointer :: f_per_bond(:, :)
  real(c_double), pointer :: wpot_per_bond(:, :, :)
            """
        s += """
  error = ERROR_NONE
  call c_f_pointer(this_cptr, this_fptr)
  call c_f_pointer(p_cptr, p)
  call c_f_pointer(nl_cptr, nl)
  if (.not. associated(this_fptr))   stop '[python_%s_energy_and_forces] *this_fptr* is NULL.'
  if (.not. associated(p))   stop '[python_%s_energy_and_forces] *p* is NULL.'
  if (.not. associated(nl))   stop '[python_%s_energy_and_forces] *nl* is NULL.'
""" % ( f90name, f90name, f90name )
        addargs = ''
        optargs = []
        if 'set_coulomb' in methods:
            addargs += 'q=q, '
        if 'mask' in features:
            s += '  if (c_associated(mask_cptr)) then\n'
            s += '     call c_f_pointer(mask_cptr, mask, [p%nat])\n'
            s += '  else\n'
            s += '     nullify(mask)\n'
            s += '  endif\n'
            optargs += ['mask']
        else:
            s += '  if (c_associated(mask_cptr)) then\n'
            s += '     RETURN_ERROR("*mask* argument present but not supported by potential %s.", error)\n' % name
            s += '  endif\n'
        if 'per_at' in features:
            s += '  if (c_associated(epot_per_at_cptr)) then\n'
            s += '     call c_f_pointer(epot_per_at_cptr, epot_per_at, [p%nat])\n'
            s += '  else\n'
            s += '     nullify(epot_per_at)\n'
            s += '  endif\n'
            s += '  if (c_associated(wpot_per_at_cptr)) then\n'
            s += '     call c_f_pointer(wpot_per_at_cptr, wpot_per_at, [3,3,p%nat])\n'
            s += '  else\n'
            s += '     nullify(wpot_per_at)\n'
            s += '  endif\n'
            optargs += ['epot_per_at', 'wpot_per_at']
        else:
            s += '  if (c_associated(epot_per_at_cptr)) then\n'
            s += '     RETURN_ERROR("*epot_per_at* argument present but not supported by potential %s.", error)\n' % name
            s += '  endif\n'
            s += '  if (c_associated(wpot_per_at_cptr)) then\n'
            s += '     RETURN_ERROR("*wpot_per_at* argument present but not supported by potential %s.", error)\n' % name
            s += '  endif\n'
        if 'per_bond' in features:
            s += '  if (c_associated(epot_per_bond_cptr)) then\n'
            s += '     call c_f_pointer(epot_per_bond_cptr, epot_per_bond, [nl%neighbors_size])\n'
            s += '  else\n'
            s += '     nullify(epot_per_bond)\n'
            s += '  endif\n'
            s += '  if (c_associated(f_per_bond_cptr)) then\n'
            s += '     call c_f_pointer(f_per_bond_cptr, f_per_bond, [3,nl%neighbors_size])\n'
            s += '  else\n'
            s += '     nullify(f_per_bond)\n'
            s += '  endif\n'
            s += '  if (c_associated(wpot_per_bond_cptr)) then\n'
            s += '     call c_f_pointer(wpot_per_bond_cptr, wpot_per_bond, [3,3,nl%neighbors_size])\n'
            s += '  else\n'
            s += '     nullify(wpot_per_bond)\n'
            s += '  endif\n'
            optargs += ['epot_per_bond', 'f_per_bond', 'wpot_per_bond']
        else:
            s += '  if (c_associated(epot_per_bond_cptr)) then\n'
            s += '     RETURN_ERROR("*epot_per_bond* argument present but not supported by potential %s.", error)\n' % name
            s += '  endif\n'
            s += '  if (c_associated(f_per_bond_cptr)) then\n'
            s += '     RETURN_ERROR("*f_per_bond* argument present but not supported by potential %s.", error)\n' % name
            s += '  endif\n'
            s += '  if (c_associated(wpot_per_bond_cptr)) then\n'
            s += '     RETURN_ERROR("*wpot_per_bond* argument present but not supported by potential %s.", error)\n' % name
            s += '  endif\n'
        if 'set_coulomb' in methods:
            s += switch_optargs('energy_and_forces_with_charges(this_fptr, p, nl, epot, f, wpot, %sierror=error)' % (addargs+'%s'), optargs)
        else:
            s += switch_optargs('energy_and_forces(this_fptr, p, nl, epot, f, wpot, %sierror=error)' % (addargs+'%s'), optargs)
        s += 'endsubroutine python_%s_energy_and_forces\n\n\n' % f90name

        f.write(s)
    
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
        s += """
void python_%s_new(void **, section_t *, section_t **);
void python_%s_free(void *);
void python_%s_register_data(void *, void *, int *);
void python_%s_init_without_parameters(void *, int *);
void python_%s_bind_to(void *, void *, void *, int *);
        """ % ( f90name, f90name, f90name, f90name, f90name )
        if 'set_coulomb' in methods:
            s += """
void python_%s_set_coulomb(void *, void *, int *);
            """ % f90name
        if 'get_dict' in methods:
            s += """
void python_%s_get_dict(void *, void *, int *);
            """ % f90name
        if 'get_per_bond_property' in methods:
            s += """
void python_%s_get_per_bond_property(void *, void *, void *, char *, double *, int *);
            """ % f90name
        s += """
void python_%s_energy_and_forces(void *, void *, void *, double *, double *, double *, int *, double *, double *, double *, double *, double *, double *, int *);
        """ % f90name

    d["prototypes"] = s

    #
    # Classes
    #

    s = "%s_class_t %s_classes[N_POTENTIAL_CLASSES] = {\n" % ( str, str )
    for f90name, f90class, name, features, methods in mods:
        s += "  {\n"
        s += "    \"%s\",\n" % name
        s += "    python_%s_new,\n" % f90name
        s += "    python_%s_free,\n" % f90name
        s += "    python_%s_register_data,\n" % f90name
        s += "    python_%s_init_without_parameters,\n" % f90name
        s += "    python_%s_bind_to,\n" % f90name
        if 'set_coulomb' in methods:
            s += "    python_%s_set_coulomb,\n" % f90name
        else:
            s += "    NULL,\n"
        if 'get_dict' in methods:
            s += "    python_%s_get_dict,\n" % f90name
        else:
            s += "    NULL,\n"
        if 'get_per_bond_property' in methods:
            s += "    python_%s_get_per_bond_property,\n" % f90name
        else:
            s += "    NULL,\n"
        s += "    python_%s_energy_and_forces,\n" % f90name
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

def write_coulomb_factory_f90(mods, str, fn):
    f = open(fn, "w")

    f.write("#include \"macros.inc\"\n\n" +
            "module %s_factory\n" % str +
            '  use libAtoms_module\n' +
            '  use particles\n' +
            '  use neighbors\n')
    for f90name, f90class, name, features, methods in mods:
        f.write('  use %s\n' % f90name)
    f.write('  implicit none\n\n' +
            'contains\n\n')

    for f90name, f90class, name, features, methods in mods:
        f.write("subroutine python_%s_new(this_cptr, cfg, m) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), intent(out) :: this_cptr\n" +
                "  type(c_ptr), value :: cfg\n" +
                "  type(c_ptr), intent(out) :: m\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  allocate(this_fptr)\n"  +
                "  if (.not. associated(this_fptr))   stop '[python_%s_new] *this_fptr* is NULL.'\n" % f90name +
                "  call register(this_fptr, cfg, m)\n" +
                "  this_cptr = c_loc(this_fptr)\n" +
                "endsubroutine python_%s_new\n\n\n" % f90name)

        f.write("subroutine python_%s_free(this_cptr) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[python_%s_free] *this_fptr* is NULL.'\n" % f90name)
        if 'del' in methods:
            f.write("  call del(this_fptr)\n")
        f.write("  deallocate(this_fptr)\n"  +
                "endsubroutine python_%s_free\n\n\n" % f90name)

        if 'register_data' in methods:
            f.write("subroutine python_%s_register_data(this_cptr, p_cptr, error) bind(C)\n" % f90name +
                    "  use, intrinsic :: iso_c_binding\n\n" +
                    "  implicit none\n\n" +
                    "  type(c_ptr), value :: this_cptr\n" +
                    "  type(c_ptr), value :: p_cptr\n" +
                    "  integer(c_int), intent(out) :: error\n\n" +
                    "  type(%s_t), pointer :: this_fptr\n" % f90name +
                    "  type(particles_t), pointer :: p\n" +
                    "  error = ERROR_NONE\n" +
                    "  call c_f_pointer(this_cptr, this_fptr)\n" +
                    "  call c_f_pointer(p_cptr, p)\n" +
                    "  if (.not. associated(this_fptr))   stop '[python_%s_register_data] *this_fptr* is NULL.'\n" % f90name +
                    "  if (.not. associated(p))   stop '[python_%s_register_data] *p* is NULL.'\n" % f90name +
                    "  call register_data(this_fptr, p, ierror=error)\n" +
                    "endsubroutine python_%s_register_data\n\n\n" % f90name)
        else:
            f.write("subroutine python_%s_register_data(this_cptr, p_cptr, error) bind(C)\n" % f90name +
                    "  use, intrinsic :: iso_c_binding\n\n" +
                    "  implicit none\n\n" +
                    "  type(c_ptr), value :: this_cptr\n" +
                    "  type(c_ptr), value :: p_cptr\n" +
                    "  integer(c_int), intent(out) :: error\n\n" +
                    "  type(%s_t), pointer :: this_fptr\n" % f90name +
                    "  type(particles_t), pointer :: p\n" +
                    "  error = ERROR_NONE\n" +
                    "  call c_f_pointer(this_cptr, this_fptr)\n" +
                    "  call c_f_pointer(p_cptr, p)\n" +
                    "  if (.not. associated(this_fptr))   stop '[python_%s_register_data] *this_fptr* is NULL.'\n" % f90name +
                    "  if (.not. associated(p))   stop '[python_%s_register_data] *p* is NULL.'\n" % f90name +
                    "endsubroutine python_%s_register_data\n\n\n" % f90name)

        f.write("subroutine python_%s_init_without_parameters(this_cptr, error) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  integer(C_INT), intent(out) :: error\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  error = ERROR_NONE\n" +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  if (.not. associated(this_fptr))   stop '[%s_init_without_parameters] *this_fptr* is NULL.'\n" % f90name)
        if 'init' in methods:
            f.write("  call init(this_fptr, error=error)\n")
        f.write("endsubroutine python_%s_init_without_parameters\n\n\n" % f90name)

        if 'set_hubbard_u' in methods:
            f.write("subroutine python_%s_set_Hubbard_U(this_cptr, p_cptr, U, error) bind(C)\n" % f90name +
                    "  use, intrinsic :: iso_c_binding\n\n" +
                    "  implicit none\n\n" +
                    "  type(c_ptr), value :: this_cptr\n" +
                    "  type(c_ptr), value :: p_cptr\n" +
                    "  real(c_double), intent(in) :: U(*)\n" +
                    "  integer(c_int), intent(out) :: error\n\n" +
                    "  type(%s_t), pointer :: this_fptr\n" % f90name +
                    "  type(particles_t), pointer :: p\n" +
                    "  error = ERROR_NONE\n" +
                    "  call c_f_pointer(this_cptr, this_fptr)\n" +
                    "  call c_f_pointer(p_cptr, p)\n" +
                    "  if (.not. associated(this_fptr))   stop '[python_%s_set_Hubbard_U] *this_fptr* is NULL.'\n" % f90name +
                    "  if (.not. associated(p))   stop '[python_%s_set_Hubbard_U] *p* is NULL.'\n" % f90name +
                    "  call set_Hubbard_U(this_fptr, p, U, error=error)\n" +
                    "endsubroutine python_%s_set_Hubbard_U\n\n\n" % f90name)

        f.write("subroutine python_%s_bind_to(this_cptr, p_cptr, nl_cptr, error) bind(C)\n" % f90name +
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
                "  if (.not. associated(this_fptr))   stop '[python_%s_bind_to] *this_fptr* is NULL.'\n" % f90name +
                "  if (.not. associated(p))   stop '[python_%s_bind_to] *p* is NULL.'\n" % f90name +
                "  if (.not. associated(nl))   stop '[python_%s_bind_to] *nl* is NULL.'\n" % f90name +
                "  call bind_to(this_fptr, p, nl, error)\n" +
                "endsubroutine python_%s_bind_to\n\n\n" % f90name)

        f.write("subroutine python_%s_energy_and_forces(this_cptr, p_cptr, nl_cptr, q, epot, f, wpot, error) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(c_ptr), value :: p_cptr\n" +
                "  type(c_ptr), value :: nl_cptr\n" +
                "  real(c_double), intent(out) :: epot\n" +
                "  real(c_double) :: q(*)\n" +
                "  real(c_double) :: f(3, *)\n" +
                "  real(c_double) :: wpot(3, 3)\n" +
                "  integer(c_int), intent(out) :: error\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  type(particles_t), pointer :: p\n" +
                "  type(neighbors_t), pointer :: nl\n" +
                "  error = ERROR_NONE\n" +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  call c_f_pointer(p_cptr, p)\n" +
                "  call c_f_pointer(nl_cptr, nl)\n" +
                "  if (.not. associated(this_fptr))   stop '[python_%s_energy_and_forces] *this_fptr* is NULL.'\n" % f90name +
                "  if (.not. associated(p))   stop '[python_%s_energy_and_forces] *p* is NULL.'\n" % f90name +
                "  if (.not. associated(nl))   stop '[python_%s_energy_and_forces] *nl* is NULL.'\n" % f90name +
                "  call energy_and_forces(this_fptr, p, nl, q, epot, f, wpot,&\n" +
                "    error)\n" +
                "endsubroutine python_%s_energy_and_forces\n\n\n" % f90name)

        f.write("subroutine python_%s_potential(this_cptr, p_cptr, nl_cptr, q, phi, error) bind(C)\n" % f90name +
                "  use, intrinsic :: iso_c_binding\n\n" +
                "  implicit none\n\n" +
                "  type(c_ptr), value :: this_cptr\n" +
                "  type(c_ptr), value :: p_cptr\n" +
                "  type(c_ptr), value :: nl_cptr\n" +
                "  real(c_double) :: q(*)\n" +
                "  real(c_double) :: phi(*)\n" +
                "  integer(c_int), intent(out) :: error\n\n" +
                "  type(%s_t), pointer :: this_fptr\n" % f90name +
                "  type(particles_t), pointer :: p\n" +
                "  type(neighbors_t), pointer :: nl\n" +
                "  error = ERROR_NONE\n" +
                "  call c_f_pointer(this_cptr, this_fptr)\n" +
                "  call c_f_pointer(p_cptr, p)\n" +
                "  call c_f_pointer(nl_cptr, nl)\n" +
                "  if (.not. associated(this_fptr))   stop '[python_%s_potential] *this_fptr* is NULL.'\n" % f90name +
                "  if (.not. associated(p))   stop '[python_%s_potential] *p* is NULL.'\n" % f90name +
                "  if (.not. associated(nl))   stop '[python_%s_potential] *nl* is NULL.'\n" % f90name +
                "  call potential(this_fptr, p, nl, q, phi, error)\n" +
                "endsubroutine python_%s_potential\n\n\n" % f90name)
    
    f.write("endmodule %s_factory\n" % str)
    f.close()

###

def write_coulomb_factory_c(mods, str, c_dispatch_template, c_dispatch_file,
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
        s += "void python_%s_new(void **, section_t *, section_t **);\n" % f90name
        s += "void python_%s_free(void *);\n" % f90name
        s += "void python_%s_register_data(void *, void *, int *);\n" % f90name
        s += "void python_%s_init_without_parameters(void *, int *);\n" % f90name
        if 'set_hubbard_u' in methods:
            s += "void python_%s_set_hubbard_u(void *, void *, double *, int *);\n" % f90name
        s += "void python_%s_bind_to(void *, void *, void *, int *);\n" % f90name
        s += "void python_%s_energy_and_forces(void *, void *, void *, double *, double *, double *, double *, int *);\n" % f90name
        s += "void python_%s_potential(void *, void *, void *, double *, double *, int *);\n" % f90name

    d["prototypes"] = s

    #
    # Classes
    #

    s = "%s_class_t %s_classes[N_COULOMB_CLASSES] = {\n" % ( str, str )
    for f90name, f90class, name, features, methods in mods:
        s += "  {\n"
        s += "    \"%s\",\n" % name
        s += "    python_%s_new,\n" % f90name
        s += "    python_%s_free,\n" % f90name
        s += "    python_%s_register_data,\n" % f90name
        s += "    python_%s_init_without_parameters,\n" % f90name
        if 'set_hubbard_u' in methods:
            s += "    python_%s_set_hubbard_u,\n" % f90name
        else:
            s += "    NULL,\n"
        s += "    python_%s_bind_to,\n" % f90name
        s += "    python_%s_energy_and_forces,\n" % f90name
        s += "    python_%s_potential,\n" % f90name
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
                    srcdir + "/c/factory.template.c", "potentials_factory_c.c",
                    srcdir + "/c/factory.template.h", "potentials_factory_c.h")
