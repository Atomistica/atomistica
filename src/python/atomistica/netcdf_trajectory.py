"""
netcdftrajectory - I/O trajectory files in the AMBER NetCDF convention

More information on the AMBER NetCDF conventions can be found at
http://ambermd.org/netcdf/. This module supports extensions to
these conventions, such as writing of additional fields and writing to
HDF5 (NetCDF-4) files.

A Python NetCDF module is required. Supported are

    netCDF4-python - http://code.google.com/p/netcdf4-python/

    scipy.io.netcdf - http://docs.scipy.org/doc/scipy/reference/io.html

    pupynere - https://bitbucket.org/robertodealmeida/pupynere/

Availability is checked in the above order of preference. Note that
scipy.io.netcdf and pupynere cannot write HDF5 NetCDF-4 files.

NetCDF files can be directly visualized using the libAtoms flavor of
AtomEye (http://www.libatoms.org/),
VMD (http://www.ks.uiuc.edu/Research/vmd/)
or Ovito (http://www.ovito.org/, starting with version 2.3).
"""

import os

import numpy as np

import ase
import ase.version

from ase.data import atomic_masses
from ase.lattice.spacegroup.cell import cellpar_to_cell, cell_to_cellpar

NC_NOT_FOUND = 0
NC_IS_NETCDF4 = 1
NC_IS_SCIPY = 2
NC_IS_PUPYNERE = 3

have_nc = NC_NOT_FOUND
# Check if we have netCDF4-python
try:
    from netCDF4 import Dataset
    have_nc = NC_IS_NETCDF4
except:
    pass

if not have_nc:
    # Check for scipy
    try:
        from scipy.io.netcdf import netcdf_file
        have_nc = NC_IS_SCIPY
    except:
        pass

if not have_nc:
    # Check for pupynere (comes with ASE)
    try:
        from ase.io.pupynere import netcdf_file
        have_nc = NC_IS_PUPYNERE
    except:
        pass


### Read/write NetCDF trajectories


class NetCDFTrajectory:
    """
    Reads/writes Atoms objects into an AMBER-stlye .nc trajectory file.
    """

    # netCDF4-python format strings to scipy.io.netcdf version numbers
    _netCDF4_to_scipy = {'NETCDF3_CLASSIC': 1, 'NETCDF3_64BIT': 2}
    _netCDF4_to_pupynere = ['NETCDF3_CLASSIC']

    # Default dimension names
    _frame_dim = 'frame'
    _spatial_dim = 'spatial'
    _atom_dim = 'atom'
    _cell_spatial_dim = 'cell_spatial'
    _cell_angular_dim = 'cell_angular'
    _label_dim = 'label'

    # Default field names. If it is a list, check for any of these names upon
    # opening. Upon writing, use the first name.
    _time_var = 'time'
    _numbers_var = ['Z', 'atom_types']
    _positions_var = 'coordinates'
    _velocities_var = 'velocities'
    _cell_origin_var = 'cell_origin'
    _cell_lengths_var = 'cell_lengths'
    _cell_angles_var = 'cell_angles'

    _default_vars = reduce(lambda x, y: x + y,
                           [_numbers_var, [_positions_var], [_velocities_var],
                           [_cell_origin_var], [_cell_lengths_var],
                           [_cell_angles_var]])

    def __init__(self, filename, mode='r', atoms=None, types_to_numbers=None,
                 double=True, netcdf_format='NETCDF3_CLASSIC'):
        """
        A NetCDFTrajectory can be created in read, write or append mode.

        Parameters:

        filename:
            The name of the parameter file.  Should end in .nc.

        mode='r':
            The mode.

            'r' is read mode, the file should already exist, and no atoms
            argument should be specified.

            'w' is write mode. The atoms argument specifies the Atoms object
            to be written to the file, if not given it must instead be given
            as an argument to the write() method.

            'a' is append mode.  It acts a write mode, except that data is
            appended to a preexisting file.

        atoms=None:
            The Atoms object to be written in write or append mode.

        types_to_numbers=None:
            Dictionary for conversion of atom types to atomic numbers when
            reading a trajectory file.

        double=True:
            Create new variable in double precision.

        netcdf_format='NETCDF3_CLASSIC':
            Format string for the underlying NetCDF file format. Only relevant
            if a new file is created. More information can be found at
            https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/File-Format.html

            'NETCDF3_CLASSIC' is the original binary format.

            'NETCDF3_64BIT' can be used to write larger files.

            'NETCDF4_CLASSIC' is HDF5 with some NetCDF limitations.

            'NETCDF4' is HDF5.
        """
        if not have_nc:
            raise RuntimeError('NetCDFTrajectory requires a NetCDF Python '
                               'module.')

        self.nc = None

        self.numbers = None
        self.pbc = None
        self.pre_observers = []   # Callback functions before write
        self.post_observers = []  # Callback functions after write
                                  # are called

        self.has_header = False
        self.set_atoms(atoms)

        self.types_to_numbers = None
        if types_to_numbers:
            self.types_to_numbers = np.array(types_to_numbers)

        # 'l' should be a valid type according to the netcdf4-python
        # documentation, but does not appear to work.
        self.dtype_conv = {'l': 'i'}
        if not double:
            self.dtype_conv.update(dict(d='f'))

        self.extra_per_frame_vars = []
        self.extra_per_file_vars = []
        # per frame atts are global quantities, not quantities stored for each
        # atom
        self.extra_per_frame_atts = []

        self.mode = mode
        self.netcdf_format = netcdf_format

        if atoms:
            self.n_atoms = len(atoms)
        else:
            self.n_atoms = None

        self._open(filename)

    def __del__(self):
        self.close()

    def _open(self, filename):
        """
        Opens the file.

        For internal use only.
        """
        if self.mode == 'a' and not os.path.exists(filename):
            self.mode = 'w'
        if have_nc == NC_IS_NETCDF4:
            self.nc = Dataset(filename, self.mode, format=self.netcdf_format)
        elif have_nc == NC_IS_SCIPY:
            if self.netcdf_format not in self._netCDF4_to_scipy:
                raise ValueError("NetCDF format '%s' not supported by "
                                 "scipy.io.netcdf." % self.netcdf_format)
            version = self._netCDF4_to_scipy[self.netcdf_format]
            if version == 1:
                # This supports older scipy.io.netcdf versions that do not
                # support the 'version' argument
                self.nc = netcdf_file(
                    filename, self.mode
                    )
            else:
                self.nc = netcdf_file(
                    filename, self.mode,
                    version=self._netCDF4_to_scipy[self.netcdf_format]
                    )
        elif have_nc == NC_IS_PUPYNERE:
            if self.netcdf_format not in self._netCDF4_to_pupynere:
                raise ValueError("NetCDF format '%s' not supported by "
                                 "ase.io.pupynere." % self.netcdf_format)
            self.nc = netcdf_file(filename, self.mode)
        else:
            # Should not happen
            raise RuntimeError('Internal error: Unknown *have_nc* value.')

        self.frame = 0
        if self.mode == 'r' or self.mode == 'a':
            self.read_header()
            self.frame = len(self)

    def set_atoms(self, atoms=None):
        """
        Associate an Atoms object with the trajectory.

        Mostly for internal use.
        """
        if atoms is not None and not hasattr(atoms, 'get_positions'):
            raise TypeError('"atoms" argument is not an Atoms object.')
        self.atoms = atoms

    def read_header(self):
        self.pbc = np.array([True, True, True])
        if not self.n_atoms:
            if have_nc == NC_IS_NETCDF4:
                self.n_atoms = len(self.nc.dimensions[self._atom_dim])
            else:
                self.n_atoms = self.nc.dimensions[self._atom_dim]
        self.numbers = self._get_variable(self._numbers_var)[:]
        if self.types_to_numbers is not None:
            self.numbers = self.types_to_numbers[self.numbers]
        self.masses = atomic_masses[self.numbers]

        for name, var in self.nc.variables.iteritems():
            # _default_vars is taken care of already
            if name not in self._default_vars:
                if len(var.dimensions) >= 2:
                    if var.dimensions[0] == self._frame_dim:
                        if var.dimensions[1] == self._atom_dim:
                            self.extra_per_frame_vars += [name]
                        else:
                            self.extra_per_frame_atts += [name]

                elif len(var.dimensions) == 1:
                    if var.dimensions[0] == self._atom_dim:
                        self.extra_per_file_vars += [name]
                    elif var.dimensions[0] == self._frame_dim:
                        self.extra_per_frame_atts += [name]

        self.has_header = True

    def write(self, atoms=None, frame=None, arrays=None, time=None):
        """
        Write the atoms to the file.

        If the atoms argument is not given, the atoms object specified
        when creating the trajectory object is used.
        """
        self._call_observers(self.pre_observers)
        if atoms is None:
            atoms = self.atoms

        if hasattr(atoms, 'interpolate'):
            # seems to be a NEB
            neb = atoms
            assert not neb.parallel
            try:
                neb.get_energies_and_forces(all=True)
            except AttributeError:
                pass
            for image in neb.images:
                self.write(image)
            return

        if not self.has_header:
            self.write_header(atoms)
        else:
            if (atoms.pbc != self.pbc).any():
                raise ValueError('Bad periodic boundary conditions!')
            elif len(atoms) != self.n_atoms:
                raise ValueError('Bad number of atoms!')
            if self.frame > 0:
                if (atoms.numbers != self.numbers).any():
                    raise ValueError('Bad atomic numbers!')
            else:
                self.numbers = atoms.get_atomic_numbers()
                self._get_variable(self._numbers_var)[:] = \
                    atoms.get_atomic_numbers()

        if frame is None:
            i = self.frame
        else:
            i = frame

        self._get_variable(self._positions_var)[i] = atoms.get_positions()
        if atoms.has('momenta'):
            self._add_velocities()
            self._get_variable(self._velocities_var)[i] = \
                atoms.get_momenta() / atoms.get_masses().reshape(-1, 1)
        a, b, c, alpha, beta, gamma = cell_to_cellpar(atoms.get_cell())
        self._get_variable(self._cell_lengths_var)[i] = [a, b, c]
        self._get_variable(self._cell_angles_var)[i] = [alpha, beta, gamma]
        if arrays is not None:
            for array in arrays:
                data = atoms.get_array(array)
                self._add_array(atoms, array, data.dtype, data.shape)
                self._get_variable(array)[i] = data
        if time is not None:
            self._get_variable(self._time_var)[i] = time

        self._call_observers(self.post_observers)
        self.frame += 1

    def write_arrays(self, atoms, frame, arrays):
        self._call_observers(self.pre_observers)
        for array in arrays:
            data = atoms.get_array(array)
            self._add_array(atoms, array, data.dtype, data.shape)
            self._get_variable(array)[frame] = data
        self._call_observers(self.post_observers)

    def _define_file_structure(self, atoms):
        if not hasattr(self.nc, 'Conventions'):
            self.nc.Conventions = 'AMBER'
        if not hasattr(self.nc, 'ConventionVersion'):
            self.nc.ConventionVersion = '1.0'
        if not hasattr(self.nc, 'program'):
            self.nc.program = 'ASE'
        if not hasattr(self.nc, 'programVersion'):
            self.nc.programVersion = ase.version.version

        if not self._frame_dim in self.nc.dimensions:
            self.nc.createDimension(self._frame_dim, None)
        if not self._spatial_dim in self.nc.dimensions:
            self.nc.createDimension(self._spatial_dim, 3)
        if not self._atom_dim in self.nc.dimensions:
            self.nc.createDimension(self._atom_dim, len(atoms))
        if not self._cell_spatial_dim in self.nc.dimensions:
            self.nc.createDimension(self._cell_spatial_dim, 3)
        if not self._cell_angular_dim in self.nc.dimensions:
            self.nc.createDimension(self._cell_angular_dim, 3)

        if not self._has_variable(self._numbers_var):
            self.nc.createVariable(self._numbers_var[0], 'i',
                                   (self._atom_dim,))
        if not self._has_variable(self._positions_var):
            self.nc.createVariable(self._positions_var, 'f4',
                                   (self._frame_dim, self._atom_dim,
                                    self._spatial_dim))
            self.nc.variables[self._positions_var].units = 'Angstrom'
            self.nc.variables[self._positions_var].scale_factor = 1.
        if not self._has_variable(self._cell_lengths_var):
            self.nc.createVariable(self._cell_lengths_var, 'd',
                                   (self._frame_dim, self._cell_spatial_dim))
            self.nc.variables[self._cell_lengths_var].units = 'Angstrom'
            self.nc.variables[self._cell_lengths_var].scale_factor = 1.
        if not self._has_variable(self._cell_angles_var):
            self.nc.createVariable(self._cell_angles_var, 'd',
                                   (self._frame_dim, self._cell_angular_dim))
            self.nc.variables[self._cell_angles_var].units = 'degree'

    def _add_velocities(self):
        if not self._has_variable(self._velocities_var):
            self.nc.createVariable(self._velocities_var, 'f4',
                                   (self._frame_dim, self._atom_dim,
                                    self._spatial_dim))
            self.nc.variables[self._positions_var].units = \
                'Angstrom/Femtosecond'
            self.nc.variables[self._positions_var].scale_factor = 1.

    def _add_array(self, atoms, array_name, type, shape):
        if not self._has_variable(array_name):
            dims = [self._frame_dim]
            for i in shape:
                if i == len(atoms):
                    dims += [self._atom_dim]
                elif i == 3:
                    dims += [self._spatial_dim]
                else:
                    raise TypeError("Don't know how to dump array of shape {0}"
                                    " into NetCDF trajectory.".format(shape))
            try:
                t = self.dtype_conv[type.char]
            except:
                t = type
            self.nc.createVariable(array_name, t, dims)

    def _get_variable(self, name):
        if isinstance(name, list):
            for n in name:
                if n in self.nc.variables:
                    return self.nc.variables[n]
            raise RuntimeError('None of the variables {0} was found in the '
                               'NetCDF trajectory.'.format(
                                   reduce(lambda x, y: x + ', ' + y, name)))
        else:
            return self.nc.variables[name]

    def _has_variable(self, name):
        if isinstance(name, list):
            for n in name:
                if n in self.nc.variables:
                    return True
            return False
        else:
            return name in self.nc.variables

    def write_header(self, atoms):
        self._define_file_structure(atoms)

        self._get_variable(self._numbers_var)[:] = \
            np.asarray(atoms.get_atomic_numbers())

    def close(self):
        """Close the trajectory file."""
        if self.nc is not None:
            self.nc.close()
            self.nc = None

    def __getitem__(self, i=-1):
        if isinstance(i, slice):
            return [self[j] for j in range(*i.indices(len(self)))]

        N = len(self)
        if 0 <= i < N:
            # Construct cell shape from cell lengths and angles
            cell = cellpar_to_cell(
                list(self.nc.variables[self._cell_lengths_var][i]) +
                list(self.nc.variables[self._cell_angles_var][i])
                )

            # Do we have a cell origin?
            if self._has_variable(self._cell_origin_var):
                origin = np.array(self.nc.variables[self._cell_origin_var][i])
            else:
                origin = np.zeros([3], dtype=float)
            origin.shape = (1, -1)

            # Compute momenta from velocities (if present)
            if self._has_variable(self._velocities_var):
                momenta = self.nc.variables[self._velocities_var][i] * \
                          self.masses.reshape(-1, 1)
            else:
                momenta = None

            # Fill info dict with additional data found in the NetCDF file
            info = {}
            for name in self.extra_per_frame_atts:
                info[name] = np.array(self.nc.variables[name][i])

            # Create atoms object
            atoms = ase.Atoms(
                positions=self.nc.variables[self._positions_var][i] - origin,
                numbers=self.numbers,
                cell=cell,
                momenta=momenta,
                masses=self.masses,
                pbc=self.pbc,
                info=info
                )

            # Attach additional arrays found in the NetCDF file
            for name in self.extra_per_frame_vars:
                atoms.set_array(name, self.nc.variables[name][i])
            for name in self.extra_per_file_vars:
                atoms.set_array(name, self.nc.variables[name][:])
            return atoms

        i = N + i
        if i < 0 or i >= N:
            raise IndexError('Trajectory index out of range.')
        return self[i]

    def __len__(self):
        if self._frame_dim in self.nc.dimensions:
            return len(self._get_variable(self._positions_var)[:])
        else:
            return 0

    def pre_write_attach(self, function, interval=1, *args, **kwargs):
        """
        Attach a function to be called before writing begins.

        function: The function or callable object to be called.

        interval: How often the function is called.  Default: every time (1).

        All other arguments are stored, and passed to the function.
        """
        if not callable(function):
            raise ValueError('Callback object must be callable.')
        self.pre_observers.append((function, interval, args, kwargs))

    def post_write_attach(self, function, interval=1, *args, **kwargs):
        """
        Attach a function to be called after writing ends.

        function: The function or callable object to be called.

        interval: How often the function is called.  Default: every time (1).

        All other arguments are stored, and passed to the function.
        """
        if not callable(function):
            raise ValueError('Callback object must be callable.')
        self.post_observers.append((function, interval, args, kwargs))

    def _call_observers(self, obs):
        """Call pre/post write observers."""
        for function, interval, args, kwargs in obs:
            if self.write_counter % interval == 0:
                function(*args, **kwargs)
