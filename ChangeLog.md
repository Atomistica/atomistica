TODO
====

- Make TabulatedEAM work for alloys. Probably useful to rename TabulatedEAM
  to FuncflEAM and TabulatedAlloyEAM to SetflEAM at this point.

Change log
==========

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
  so far independent of ASEs PBS setting.)

v0.2.0
------

- Added support for non-orthogonal tight-binding (NOTB) with charge
  self-consistency (SCC).

- Added Coulomb solvers (required for SCC-NOTB).

v0.1.3
------

- Compatibility with IBM XL compiler suite on the BlueGene architecture.

