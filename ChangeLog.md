Change log
==========

v1.0.4 (24Apr22)
----------------

- Automatically use correct GCC switches for legacy Fortran code

v1.0.3 (06Aug21)
----------------

- Bug fix: SlaterCharges and GaussianCharges could case segmentation fault
  because of uninitialized value

v1.0.2 (19Jul21)
----------------

- Fixed regression in makefile.inc

v1.0.1 (19Jul21)
----------------

- Minor fixes for the Intel compiler suite
- Docker and Singularity containers

v1.0.0 (25Feb21)
-----------------

- Drop support for Python < 3.6 (#34)
- Bug fix: Reading EAM tables
- Bug fix: Python bindings should recalculate only when necessary 

v0.10.2 (28May20)
-----------------
- Update Python interface for recent changes in ASE
- Bug fix: Wrong neighbor list construction for small negative positions (#30)
- Bug fix: Default system of units in DFT-D3 when using Python or LAMMPS (#27)

v0.10.1 (18Nov19)
-----------------

- Select elements in SlidingT

v0.10.0 (17Nov19)
-----------------

- DFT-D3 dispersion correction
- Variable-charge model

v0.9.2 (5Nov19)
---------------

- OpenMP parallelization of the grid portion of PME

v0.9.1 (4Nov19)
---------------

- OpenMP parallelization of LJCut

v0.9.0 (3Nov19)
---------------

- Variable charge model

v0.8.0 (14Apr19)
----------------

- Support for DFTB3.

v0.7.0 (9Dec18)
---------------

- Updated to support latest ASE master.
- Fixed random hangs when running with LAMMPS.
- Bug fix in C-H interaction of screened REBO2 potential when running
  within LAMMPS.
  
v0.6.0 (18May18)
----------------

- Updated for ASE 3.15.
- Bug fix in SCC-DFTB.

v0.5.6 (10May18)
----------------

- Fixed memory leak.

v0.5.5 (10Apr18)
----------------

- Improved introspection of the electronic structure obtained in tight-binding
  calculations.

v0.5.4 (6Oct17)
---------------

- Fixed segfault in tight-binding materials database on Mac OS X.

v0.5.3 (4Jul17)
---------------

- Support for ASE's Atoms.celldisp parameter
- Bug fix: Atoms outside the simulation cell are now treated correctly in
  periodic and nonperiodic cells.

v0.5.2 (12Jun17)
----------------

- Removed LAPACK dependency for everything that does not use tight-binding
  (in particular the LAMMPS moduls)
- Updated LAMMPS pair_style for LAMMPS >= 07Sep16
- Bug fix: Proper inclusion of numpy extra_link_args in setup.py
- Bug fix: Fixed problem with tilted orthorhombic cells

v0.5.1 (28Feb17)
----------------

- Bug fix: Set element charge from skf file. Only affects tight-binding runs
  with Slater-Koster tables from dftb.org

v0.5.0 (27Feb17)
----------------

- Particle mesh Ewald Coulomb solver
- Regression fix in empirical bond-order potentials (dimers were handle
  incorrectly) - regression introduced in 0.4.4

v0.4.6 (11Nov16)
----------------

- More robust EAM (low and high density configurations)
- More robust version extraction for LAMMPS and standalone code (uses
  versioneer)

v0.4.5 (19Aug16)
----------------

- Maintenance: Fixed problem with passing c_null_ptr through c_f_pointer when
  using XLF

v0.4.4 (19Jul16)
----------------

- Bug fix: Occasional NaNs in bond-order potentials
- Maintenance: Fixed segfault on BlueGene

v0.4.3 (20Jun16)
----------------

- Preparation for PyPI release
- Minor fixes to tight-binding solver (when used with MPI domain
  decomposition)

v0.4.2 (24Mar16)
----------------

- Support for new ASE calculator interface (thanks to James Kermode)

v0.4.1 (18Oct15)
----------------

- Regression fix: Charge extrapolation did not work because charges were
  overriden with initial charges at every step. This was introduced when
  changing to the new ASE Calculator class.

v0.4.0 (12Oct15)
----------------

- Python 3 compatibility.

- Python interface based on new ASE Calculator base class.

- Tight binding: Removed necessity for an 'elements.dat' file.
  Code auto-detects elements from files found in the database directory.

- Tight binding: Added support for d, sd and pd electronic configurations.

- LAMMPS interface automatically checks git fingerprint to ensure
  compatibility between pair_style and Atomistica library.

- Fixed proper stopping of timer when errors are raised or passed in some
  parts of the code.

v0.3.2 (13Apr15)
----------------

- General mechanism for object introspection from Python

- Exposed NOTB internals to Python (Hamiltonian, overlap and density matrices,
  eigenvalues and eigenvectors if LAPACK solver is used)

- Made NOTB per-bond analysis tools available from Python (added
  get_per_bond_property method)

- Some bug fixes to parts of the standalone code

v0.3.1 (28Feb15)
----------------

- Implemented charge extrapolation for SCC NOTB

- Fixed a couple of OpenMP related regression that broke compilation.

- Fixed a regression that lead to wrong unit conversion for Hubbard-Us.

- Fixed a regression that lead to OutputEnergy not being called.

- Coulomb solvers now have an energy_and_forces rather than a
  potential_and_field function.

v0.3.0 (18Feb15)
----------------

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

v0.2.5 (16Dec14)
----------------

- Corrected another bug regarding handling of PBCs in ASE.

v0.2.4 (6Nov14)
---------------

- Corrected handling of PBCs in ASE.

v0.2.3 (26Sep14)
----------------

- Fixed buffer overflow bug in atom type handling in LAMMPS interface.

v0.2.2 (13Aug14)
----------------

- Added harmonic spring (with cutoff) potential.
- Added double harmonic potential. (Two harmonic springs at different distance,
  useful to model an SC solid that is isotropically elastic.)

v0.2.1 (18Jul14)
----------------

- Added proper handling of PBCs to ASE interface. (All systems were periodic
  so far independent of ASEs PBC setting.)

v0.2.0 (24Jun14)
----------------

- Added support for non-orthogonal tight-binding (NOTB) with charge
  self-consistency (SCC).

- Added Coulomb solvers (required for SCC-NOTB).

v0.1.3 (7Jan14)
---------------

- Compatibility with IBM XL compiler suite on the BlueGene architecture.
