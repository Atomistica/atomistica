import numpy as np

class JoinCalculators:
    """Create a joint calculator, i.e. one can augment a DFT calculation
       by adding a classical potential for van-der-Waals interactions.

       Potential energies, forces, etc. are simply summed up.
    """

    def __init__(self, calcs):
        self.calcs = calcs


    def get_forces(self, a):
        """Calculate atomic forces."""
        f = np.zeros( [ len(a), 3 ], dtype=float )
        for c in self.calcs:
            f += c.get_forces(a)
        return f


    def get_magnetic_moments(self, a):
        """Get calculated local magnetic moments."""
        raise NotImplementedError


    def get_potential_energy(self, a):
        """Calculate potential energy."""
        e = 0.0
        for c in self.calcs:
            e += c.get_potential_energy(a)
        return e


    def get_potential_energies(self, a):
        """Calculate the potential energies of all the atoms."""
        raise NotImplementedError


    def get_spin_polarized(self):
        """Get calculated total magnetic moment."""
        raise NotImplementedError


    def get_stress(self, a):
        """Calculate stress tensor."""
        s = np.zeros( 6, dtype=float )
        for c in self.calcs:
            s += c.get_stress(a)
        return s


    def get_stresses(self, a):
        """Calculate the stress-tensor of all the atoms."""
        raise NotImplementedError


    def set_atoms(self, a):
        """Assign an atoms object."""
        for c in self.calcs:
            if hasattr(c, "set_atoms"):
                c.set_atoms(a)

