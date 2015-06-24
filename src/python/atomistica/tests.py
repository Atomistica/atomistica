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

"""
Potential test suite.
"""

from __future__ import print_function

from math import sqrt

import numpy as np

import ase
import ase.constraints
from ase.optimize import FIRE, QuasiNewton
from ase.units import GPa
import ase.lattice.cubic as cubic

###

Jm2 = 1e23/ase.units.kJ

###

def test_forces(atoms, dx=1e-6):
    """Compute forces and compare to forces computed numerically from a 
       finite differences approach.
    """

    f0  = atoms.get_forces().copy()
    ffd = f0.copy()

    for a in atoms:
        r0  = a.position.copy()

        a.x = r0[0]-dx
        ex1 = atoms.get_potential_energy()
        a.x = r0[0]+dx
        ex2 = atoms.get_potential_energy()
        a.x = r0[0]

        a.y = r0[1]-dx
        ey1 = atoms.get_potential_energy()
        a.y = r0[1]+dx
        ey2 = atoms.get_potential_energy()
        a.y = r0[1]

        a.z = r0[2]-dx
        ez1 = atoms.get_potential_energy()
        a.z = r0[2]+dx
        ez2 = atoms.get_potential_energy()
        a.z = r0[2]

        ffd[a.index, 0] = -(ex2-ex1)/(2*dx)
        ffd[a.index, 1] = -(ey2-ey1)/(2*dx)
        ffd[a.index, 2] = -(ez2-ez1)/(2*dx)

    df    = ffd-f0
    absdf = np.sum(df*df, axis=1)

    return ffd, f0, np.max(absdf)


def test_virial(atoms, de=1e-6):
    """Compute virial and compare to virial computed numerically from a 
       finite differences approach.
    """

    s0  = atoms.get_stress().copy()
    V0  = atoms.get_volume()
    sfd = np.zeros([ 3, 3 ])
    c0  = atoms.get_cell().copy()

    un      = np.zeros([3,3])
    un[0,0] = 1.0
    un[1,1] = 1.0
    un[2,2] = 1.0


    for i in range(3):
        for j in range(3):
            c         = c0.copy()
            eps       = un.copy()

            eps[i, j] = un[i, j]-de
            c         = np.dot(c0, eps)
            atoms.set_cell(c, scale_atoms=True)
            e1  = atoms.get_potential_energy()

            eps[i, j] = un[i, j]+de
            c         = np.dot(c0, eps)
            atoms.set_cell(c, scale_atoms=True)
            e2  = atoms.get_potential_energy()

            sfd[i, j] = (e2-e1)/(2*de)

    sfd = np.array( [ sfd[0,0], sfd[1,1], sfd[2,2], (sfd[1,2]+sfd[2,1])/2,
                      (sfd[0,2]+sfd[2,0])/2, (sfd[0,1]+sfd[1,0])/2 ] )/V0

    return sfd, s0, np.max(sfd-s0)


def test_potential(atoms, dq=1e-6):
    """
    Compute electrostatic potential and compare to potential computed
    numerically from a finite differences approach.
    """

    p0  = atoms.get_calculator().get_electrostatic_potential().copy()
    pfd = p0.copy()

    for a in atoms:
        q0       = a.charge

        a.charge = q0-dq
        eq1      = atoms.get_potential_energy()
        a.charge = q0+dq
        eq2      = atoms.get_potential_energy()
        a.charge = q0

        pfd[a.index] = (eq2-eq1)/(2*dq)

    dp    = pfd-p0
    absdp = np.sum(dp*dp)

    return pfd, p0, np.max(absdp)


def cubic_elastic_constants(a, Minimizer=None, fmax=0.025, eps=0.001):
    r0  = a.get_positions().copy()

    cell  = a.get_cell()
    sxx0, syy0, szz0, syz0, szx0, sxy0  = a.get_stress()

    ## C11
    T = np.diag( [ eps, 0.0, 0.0 ] )
    a.set_cell( np.dot(np.eye(3)+T, cell.T).T, scale_atoms=True )
    if Minimizer is not None:
        Minimizer(a, logfile=None).run(fmax=fmax)
    sxx11, syy11, szz11, syz11, szx11, sxy11  = a.get_stress()

    C11  = (sxx11-sxx0)/eps

    ## C12 (C)
    T = np.diag( [ eps, -eps/2, -eps/2 ] )
    a.set_cell( np.dot(np.eye(3)+T, cell.T).T, scale_atoms=True )
    if Minimizer is not None:
        Minimizer(a, logfile=None).run(fmax=fmax)
    sxx12, syy12, szz12, syz12, szx12, sxy12  = a.get_stress()

    Cp   = ((sxx12-sxx0)-(syy12-syy0))/(3*eps)
    C12  = C11-2*Cp

    ## C44
    T = np.array( [ [ 0.0, 0.5*eps, 0.5*eps ], [ 0.5*eps, 0.0, 0.5*eps ], [ 0.5*eps, 0.5*eps, 0.0 ] ] )
    a.set_cell( np.dot(np.eye(3)+T, cell.T).T, scale_atoms=True )
    if Minimizer is not None:
        Minimizer(a, logfile=None).run(fmax=fmax)
    sxx44, syy44, szz44, syz44, szx44, sxy44  = a.get_stress()

    C44  = (syz44+szx44+sxy44-syz0-szx0-sxy0)/(3*eps)

    a.set_cell( cell, scale_atoms=True )
    a.set_positions(r0)

    B  = (C11+2*C12)/3

    return ( C11, C12, C44, B, Cp )


def orthorhombic_elastic_constants(a, Minimizer=None, fmax=0.025, eps=0.001):
    if Minimizer is not None:
        Minimizer(a, logfile='min.log').run(fmax=fmax)

    r0    = a.get_positions().copy()

    cell  = a.get_cell()
    s0    = a.get_stress()

    ## C11
    C11  = [ ]
    for i in range(3):
        a.set_cell(cell, scale_atoms=True)
        a.set_positions(r0)
        
        T        = np.zeros( (3,3) )
        T[i, i]  = eps
        a.set_cell( np.dot(np.eye(3)+T, cell), scale_atoms=True )
        if Minimizer is not None:
            Minimizer(a, logfile='min.log').run(fmax=fmax)
        s  = a.get_stress()
            
        C11 += [ (s[i]-s0[i])/eps ]

    ## C12 (C)
    Cp   = [ ] 
    C12  = [ ]
    for i in range(3):
        a.set_cell(cell, scale_atoms=True)
        a.set_positions(r0)
        
        T       = np.zeros( (3, 3) )
        j       = (i+1)%3
        k       = (i+2)%3
        T[j,j]  = eps
        T[k,k]  = -eps
        a.set_cell( np.dot(np.eye(3)+T, cell), scale_atoms=True )
        if Minimizer is not None:
            Minimizer(a, logfile='min.log').run(fmax=fmax)
        s  = a.get_stress()

        Cp  += [ ((s[j]-s0[j])-(s[k]-s0[k]))/(4*eps) ]

    ## C44
    C44  = [ ]
    for i in range(3):
        a.set_cell(cell, scale_atoms=True)
        a.set_positions(r0)

        T        = np.zeros( (3, 3) )
        j        = (i+1)%3
        k        = (i+2)%3
        T[j, k]  = eps
        T[k, j]  = eps
        a.set_cell( np.dot(np.eye(3)+T, cell), scale_atoms=True )
        if Minimizer is not None:
            Minimizer(a, logfile='min.log').run(fmax=fmax)
        s  = a.get_stress()
        #sxx44, syy44, szz44, syz44, szx44, sxy44  = a.get_stress()

        #C44  = (syz44+szx44+sxy44-syz0-szx0-sxy0)/(3*eps)
        C44 += [ (s[3+i]-s0[3+i])/(2*eps) ]

    a.set_cell( cell, scale_atoms=True )
    a.set_positions(r0)

    C11  = np.array(C11)
    Cp   = np.array(Cp)
    C44  = np.array(C44)
    
    C12  = C11-2*Cp

    return ( C11, C12, C44, Cp )


def test_cubic_elastic_constants(mats, pot, par=None, sx=1, dev_thres=5,
                                 test=None):
    nok = 0
    nfail = 0
    try:
        potname = pot.__name__
    except:
        potname = pot.__class__.__name__
    if test is None:
        print('--- %s ---' % potname)
    if par is not None:
        if test is None and '__ref__' in par:
            print('    %s' % par['__ref__'])
        c  = pot(**par)
    else:
        c  = pot
    for imat in mats:
        t_Ec = t_a0 = t_C11 = t_C12 = t_C44 = t_C440 = t_B = t_Cp = None
        if isinstance(imat, tuple):
            name, a, t_Ec, t_a0, t_C11, t_C12, t_C44, t_B, t_Cp = imat
        else:
            name = imat['name']
            a = imat['struct']
            try:
                t_Ec = float(imat['Ec'])
            except:
                t_Ec = None
            try:
                t_a0 = float(imat['a0'])
            except:
                t_a0 = None
            try:
                t_C11 = float(imat['C11'])
            except:
                t_C11 = None
            try:
                t_C12 = float(imat['C12'])
            except:
                t_C12 = None
            try:
                t_C44 = float(imat['C44'])
            except:
                t_C44 = None
            try:
                t_C440 = float(imat['C440'])
            except:
                t_C440 = None
            try:
                t_B = float(imat['B'])
            except:
                t_B = None
            try:
                t_Cp = float(imat['Cp'])
            except:
                t_Cp = None

        errmsg = 'potential: %s; material: %s' % (potname, name)

        a.translate([0.1, 0.1, 0.1])
        a.set_scaled_positions(a.get_scaled_positions())
        a.set_calculator(c)

        FIRE(
            ase.constraints.StrainFilter(a, mask=[1,1,1,0,0,0]),
            logfile=None).run(fmax=0.0001)

        #ase.io.write('%s.cfg' % name, a)

        #
        # Ec
        #

        Ec = a.get_potential_energy()/len(a)
        if t_Ec is None:
            if test is None:
                print('%10s: Ec   = %10.3f eV' % ( name, Ec ))
        else:
            t_Ec = float(t_Ec)
            dev = (Ec + t_Ec)*100/t_Ec
            if test is None:
                print('%10s: Ec   = %10.3f eV  (%10.3f eV - %7.2f %%)' % \
                    ( name, Ec, t_Ec, dev ))
            if test is None:
                if abs(dev) > dev_thres:
                    print('            --- Warning: Property off by more than '\
                        '%i %%.' % dev_thres)
                    nfail += 1
                else:
                    nok += 1
            else:
                test.assertTrue(abs(dev) < dev_thres, msg=errmsg)

        #
        # a0
        #

        c1, c2, c3  = a.get_cell()
        a0          = sqrt(np.dot(c1, c1))/sx
        if t_a0 is None:
            if test is None:
                print('            a0   = %10.3f A ' % a0)
        else:
            t_a0 = float(t_a0)
            dev = (a0 - t_a0)*100/t_a0
            if test is None:
                print('            a0   = %10.3f A   (%10.3f A - %7.2f %%)' % \
                    ( a0, t_a0, dev ))
                if abs(dev) > dev_thres:
                    print('            --- Warning: Property off by more than '\
                        '%i %%.' % dev_thres)
                    nfail += 1
                else:
                    nok += 1
            else:
                test.assertTrue(abs(dev) < dev_thres, msg=errmsg)

        C11, C12, C44, B, Cp       = cubic_elastic_constants(a, eps=1e-6)
        C11r, C12r, C44r, Br, Cpr  = cubic_elastic_constants(
            a, Minimizer=QuasiNewton, fmax=1e-8, eps=0.001)

        #
        # C11
        #

        if t_C11 is None:
            if test is None:
                print('            C11  = %10.4f GPa' % (C11/GPa))
        else:
            t_C11 = float(t_C11)
            dev = (C11/GPa - t_C11)*100/t_C11
            if test is None:
                print('            C11  = %10.4f GPa (%10.4f GPa - ' \
                      '%7.2f%%)' % (C11/GPa, t_C11, dev))
                if abs(dev) > dev_thres:
                    print('            --- Warning: Property off by more than '\
                        '%f %%.' % dev_thres)
                    nfail += 1
                else:
                    nok += 1
            else:
                test.assertTrue(abs(dev) < dev_thres, msg=errmsg)

        #
        # C12
        #

        if t_C12 is None:
            if test is None:
                print('            C12  = %10.4f GPa GPa' %  (C12/GPa))
        else:
            t_C12 = float(t_C12)
            dev = (C12/GPa - t_C12)*100/t_C12
            if test is None:
                print('            C12  = %10.4f GPa (%10.4f GPa ' \
                    '- %7.2f %%)' % (C12/GPa, t_C12, dev))
                if abs(dev) > dev_thres:
                    print('            --- Warning: Property off by more than '\
                        '%f %%.' % (dev_thres))
                    nfail += 1
                else:
                    nok += 1
            else:
                test.assertTrue(abs(dev) < dev_thres, msg=errmsg)

        #
        # C44
        #

        if t_C44 is None:
            if test is None:
                print('            C44  = %10.4f GPa' % (C44r/GPa))
        else:
            t_C44 = float(t_C44)
            dev = (C44r/GPa - t_C44)*100/t_C44
            if test is None:
                print('            C44  = %10.4f GPa (%10.4f GPa - ' \
                    '%7.2f %%)' % ( C44r/GPa, t_C44, dev ))
            if test is None:
                if abs(dev) > dev_thres:
                    print('            --- Warning: Property off by more than '\
                        '%f %%.' % (dev_thres))
                    nfail += 1
                else:
                    nok += 1
            else:
                test.assertTrue(abs(dev) < dev_thres, msg=errmsg)

        #
        # C440
        #

        if t_C440 is None:
            if test is None:
                print('            C440 = %10.4f GPa' % (C44/GPa))
        else:
            t_C440 = float(t_C440)
            dev = (C44/GPa - t_C440)*100/t_C440
            if test is None:
                print('            C440 = %10.4f GPa (%10.4f GPa - ' \
                      '%7.2f %%)' % (C44/GPa, t_C440, dev ))
            if test is None:
                if abs(dev) > dev_thres:
                    print('            --- Warning: Property off by more than '\
                        '%f %%.' % (dev_thres))
                    nfail += 1
                else:
                    nok += 1
            else:
                test.assertTrue(abs(dev) < dev_thres, msg=errmsg)

        #
        # B
        #

        if t_B is None:
            if test is None:
                print('            B    = %10.4f GPa' % (B/GPa))
        else:
            t_B = float(t_B)
            dev = (B/GPa - t_B)*100/t_B
            if test is None:
                print('            B    = %10.4f GPa (%10.4f GPa ' \
                    '- %7.2f %%)' % (B/GPa, t_B, dev))
            if test is None:
                if abs(dev) > dev_thres:
                    print('            --- Warning: Property off by more than '\
                        '%f %%.' % (dev_thres))
                    nfail += 1
                else:
                    nok += 1
            else:
                test.assertTrue(abs(dev) < dev_thres, msg=errmsg)

        #
        # Cp
        #

        if t_Cp is None:
            if test is None:
                print('            Cp   = %10.4f GPa' % (Cp/GPa))
        else:
            t_Cp = float(t_Cp)
            dev = (Cp/GPa - t_Cp)*100/t_Cp
            if test is None:
                print('            Cp   = %10.4f GPa (%10.4f GPa ' \
                    '- %7.2f %%)'% (Cp/GPa, t_Cp, dev))
            if test is None:
                if abs(dev) > dev_thres:
                    print('            --- Warning: Property off by more than '\
                        '%f %%.' % (dev_thres))
                    nfail += 1
                else:
                    nok += 1
            else:
                test.assertTrue(abs(dev) < dev_thres, msg=errmsg)

    return nok, nfail


def test_hexagonal_elastic_constants(mats, pot, par=None, sx=1, dev_thres=5,
                                     test=None):
    try:
        potname = pot.__name__
    except:
        potname = pot.__class__.__name__
    if test is None:
        print('--- %s ---' % potname)
    if par is not None:
        if test is None and '__ref__' in par:
            print('    %s' % par['__ref__'])
        c  = pot(**par)
    else:
        c  = pot
    for imat in mats:
        if isinstance(imat, tuple):
            name, a, t_Ec, t_a0, t_c0
        else:
            name = imat['name']
            a = imat['struct']
            try:
                t_Ec = float(imat['Ec'])
            except:
                t_Ec = None
            try:
                t_a0 = float(imat['a0'])
            except:
                t_a0 = None
            try:
                t_c0 = float(imat['c0'])
            except:
                t_c0 = None
        a.translate([0.1, 0.1, 0.1])
        a.set_scaled_positions(a.get_scaled_positions()%1.0)
        a.set_calculator(c)

        FIRE(
            ase.constraints.StrainFilter(a, mask=[1,1,0,0,0,0]),
            logfile=None).run(fmax=0.0001)

        ase.io.write('%s.cfg' % name, a)

        Ec = a.get_potential_energy()/len(a)
        if t_Ec is None:
            print('%10s: Ec  = %10.3f eV' % ( name, Ec ))
        else:
            dev = (Ec + t_Ec)*100/t_Ec
            print('%10s: Ec  = %10.3f eV  (%10.3f eV - %7.2f %%)' % ( name, Ec, t_Ec, dev ))
            if abs(dev) > dev_thres:
                print('            --- Warning: Property off by more than %i %%.' % dev_thres)

        c1, c2, c3  = a.get_cell()
        a0          = sqrt(np.dot(c1, c1))/sx
        b0          = sqrt(np.dot(c2, c2))/sx
        c0          = sqrt(np.dot(c3, c3))/sx
        a0 = (a0/sqrt(3.0)+b0)/2
        #a0 /= sqrt(3.0)
        if t_a0 is None:
            print('            a0  = %10.3f A ' % a0)
        else:
            dev = (a0 - t_a0)*100/t_a0
            print('            a0  = %10.3f A   (%10.3f A - %7.2f %%)' % ( a0, t_a0, dev ))
            if abs(dev) > dev_thres:
                print('            --- Warning: Property off by more than %i %%.' % dev_thres)
        if t_c0 is None:
            print('            c0  = %10.3f A ' % c0)
        else:
            dev = (c0 - t_c0)*100/t_c0
            print('            c0  = %10.3f A   (%10.3f A - %7.2f %%)' % ( c0, t_c0, dev ))
            if abs(dev) > dev_thres:
                print('            --- Warning: Property off by more than %i %%.' % dev_thres)


def test_surface_energies(mats, pot, par=None, sx=1, vacuum=10.0, find_a0=True,
                          dev_thres=5, test=None, dump=False):
    try:
        potname = pot.__name__
    except:
        potname = pot.__class__.__name__
    if test is None:
        print('--- %s ---' % potname)
    if par is not None:
        if test is None and '__ref__' in par:
            print('    %s' % par['__ref__'])
        c = pot(**par)
    else:
        c = pot
    for imat in mats:
        t_Es_u = t_Es_r = t_Es_u_Jm2 = t_Es_r_Jm2 = None
        if isinstance(imat, tuple):
            name, a = imat
        else:
            name = imat['name']
            a = imat['struct']
            try:
                t_Es_u = float(imat['u'])
            except:
                t_Es_u = None
            try:
                t_Es_r = float(imat['r'])
            except:
                t_Es_r = None
            try:
                t_Es_u_Jm2 = float(imat['u_Jm2'])
            except:
                t_Es_u_Jm2 = None
            try:
                t_Es_r_Jm2 = float(imat['r_Jm2'])
            except:
                t_Es_r_Jm2 = None

        errmsg = 'potential: %s; material: %s' % (potname, name)

        bulk = None
        if type(a) == tuple:
            bulk, a = a
            bulk.translate([0.1, 0.1, 0.1])
            bulk.set_scaled_positions(bulk.get_scaled_positions())
            bulk.set_calculator(c)

        a.translate([0.1, 0.1, 0.1])
        a.set_scaled_positions(a.get_scaled_positions())
        a.set_calculator(c)

        if bulk is None:
            bulk = a

        if find_a0:
            FIRE(
                ase.constraints.StrainFilter(bulk, mask=[1,1,1,0,0,0]),
                logfile=None).run(fmax=0.0001)

        Ebulk = bulk.get_potential_energy()
        if test is None:
            print('%-20s: Ec            = %10.3f eV' % (name, Ebulk/len(a)))

        cx, cy, cz = bulk.get_cell().diagonal()
        a.set_cell([cx,cy,cz], scale_atoms=True)
        a.set_cell([cx,cy,cz+vacuum])

        Eunrelaxed = a.get_potential_energy()
        # Factor of two because there are two surfaces!
        Es = ( Eunrelaxed - Ebulk ) / 2
        Es_Jm2 = Es*Jm2/(cx*cy)
        Es /= sx*sx
        if test is None:
            print('                      Es,unrelaxed  = %10.3f eV/cell ' \
                  '(%10.3f J/m^2)' % (Es, Es_Jm2))
        else:
            if t_Es_u is not None:
                dev = (Es - t_Es_u)*100/t_Es_u
                test.assertTrue(abs(dev) < dev_thres,
                                msg='Es,unrelaxed; '+errmsg)
            if t_Es_u_Jm2 is not None:
                dev = (Es_Jm2 - t_Es_u_Jm2)*100/t_Es_u_Jm2
                test.assertTrue(abs(dev) < dev_thres,
                                msg='Es,relaxed; '+errmsg)

        FIRE(a, logfile=None).run(fmax=0.005)

        Erelaxed = a.get_potential_energy()
        # Factor of two because there are two surfaces!
        Es = ( Erelaxed - Ebulk ) / 2
        Es_Jm2 = Es*Jm2/(cx*cy)
        Es /= sx*sx
        if test is None:
            print('                      Es,relaxed    = %10.3f eV/cell ' \
                  '(%10.3f J/m^2)' % (Es, Es_Jm2))
        else:
            if t_Es_r is not None:
                dev = (Es - t_Es_r)*100/t_Es_r
                test.assertTrue(abs(dev) < dev_thres,
                                msg='Es,unrelaxed (J/m^2); '+errmsg)
            if t_Es_r_Jm2 is not None:
                dev = (Es_Jm2 - t_Es_r_Jm2)*100/t_Es_r_Jm2
                test.assertTrue(abs(dev) < dev_thres, 
                                msg='Es,relaxed (J/m^2); '+errmsg)

        if dump:
            ase.io.write('%s-%s.cfg' % ( potname, name ), a)
