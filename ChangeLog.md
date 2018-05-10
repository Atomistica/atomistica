TODO
====

- Make TabulatedEAM work for alloys. Probably useful to rename TabulatedEAM
  to FuncflEAM and TabulatedAlloyEAM to SetflEAM to distinguish by file
  format.

Change log
==========

v0.5.6
------

- Fixed memory leak.

v0.5.5
------

- Improved introspection of the electronic structure obtained in tight-binding
  calculations.

v0.5.4
------

- Fixed segfault in tight-binding materials database on Mac OS X.

v0.5.3
------

- Support for ASE's Atoms.celldisp parameter
- Bug fix: Atoms outside the simulation cell are now treated correctly in
  periodic and nonperiodic cells.

v0.5.2
------

- Removed LAPACK dependency for everything that does not use tight-binding
  (in particular the LAMMPS moduls)
- Updated LAMMPS pair_style for LAMMPS >= 07Sep16
- Bug fix: Proper inclusion of numpy extra_link_args in setup.py
- Bug fix: Fixed problem with tilted orthorhombic cells

v0.5.1
------

- Bug fix: Set element charge from skf file. Only affects tight-binding runs
  with Slater-Koster tables from dftb.org

v0.5.0
------

- Particle mesh Ewald Coulomb solver
- Regression fix in empirical bond-order potentials (dimers were handle
  incorrectly) - regression introduced in 0.4.4

v0.4.6
------

- More robust EAM (low and high density configurations)
- More robust version extraction for LAMMPS and standalone code (uses
  versioneer)

v0.4.5
------

- Maintenance: Fixed problem with passing c_null_ptr through c_f_pointer when
  using XLF

v0.4.4
------

- Bug fix: Occasional NaNs in bond-order potentials
- Maintenance: Fixed segfault on BlueGene

v0.4.3
------

- Preparation for PyPI release
- Minor fixes to tight-binding solver (when used with MPI domain
  decomposition)

v0.4.2
------

- Support for new ASE calculator interface (thanks to James Kermode)

v0.4.1
------

- Regression fix: Charge extrapolation did not work because charges were
  overriden with initial charges at every step. This was introduced when
  changing to the new ASE Calculator class.

v0.4.0
------

- Python 3 compatibility.

- Python interface based on new ASE Calculator base class.

- Tight binding: Removed necessity for an 'elements.dat' file.
  Code auto-detects elements from files found in the database directory.

- Tight binding: Added support for d, sd and pd electronic configurations.

- LAMMPS interface automatically checks git fingerprint to ensure
  compatibility between pair_style and Atomistica library.

- Fixed proper stopping of timer when errors are raised or passed in some
  parts of the code.

v0.3.2
------

- General mechanism for object introspection from Python

- Exposed NOTB internals to Python (Hamiltonian, overlap and density matrices,
  eigenvalues and eigenvectors if LAPACK solver is used)

- Made NOTB per-bond analysis tools available from Python (added
  get_per_bond_property method)

- Some bug fixes to parts of the standalone code

v0.3.1
------

- Implemented charge extrapolation for SCC NOTB

- Fixed a couple of OpenMP related regression that broke compilation.

- Fixed a regression that lead to wrong unit conversion for Hubbard-Us.

- Fixed a regression that lead to OutputEnergy not being called.

- Coulomb solvers now have an energy_and_forces rather than a
  potential_and_field function.

v0.3.0
------

- Added standalone molecular dynamics code. Source is in src/standalone and
  Makefiles are in build_standalone.

- Dipatch generator now parses the additional 'features' key in the @meta
  data. Features can be 'mask', 'per_at' and 'per_bond'. This enables passing
  of *mask*, *epot/wpot_per_at* and *epot/wpot_per_bond* parameters to the
  potential. (These parameters can be omitted from the interface if unused now.)

- Added *mask* parameter that can be used to turn on/off contribution from
  individual atoms. Contributions are additive, i.e.
  epot(mask) + epot(.not. mask) gives the total energy.

- Fixed TabulatedEAM. Never worked properly but does now for unary systems.

- Size of internal neighbor lists (EAM potentials and BOPs) is now computed
  on-the-fly.

v0.2.5
------

- Corrected another bug regarding handling of PBCs in ASE.

v0.2.4
------

- Corrected handling of PBCs in ASE.

v0.2.3
------

- Fixed buffer overflow bug in atom type handling in LAMMPS interface.

v0.2.2
------

- Added harmonic spring (with cutoff) potential.
- Added double harmonic potential. (Two harmonic springs at different distance,
  useful to model an SC solid that is isotropically elastic.)

v0.2.1
------

- Added proper handling of PBCs to ASE interface. (All systems were periodic
  so far independent of ASEs PBC setting.)

v0.2.0
------

- Added support for non-orthogonal tight-binding (NOTB) with charge
  self-consistency (SCC).

- Added Coulomb solvers (required for SCC-NOTB).

v0.1.3
------

- Compatibility with IBM XL compiler suite on the BlueGene architecture.
