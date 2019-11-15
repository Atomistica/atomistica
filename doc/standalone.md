Running the standalone code
===========================

Introduction
------------

The standalone code is called MDCORE. It is able to run molecular dynamics
calculations or perform optimizations of a single structure. Unlike other
codes like LAMMPS, each call to executable will run either dynamics or
optimization under a fixed set of conditions. For example, running an
optimization followed by an NVT equilibration and an NVE production run
requires three calls of the executable.

Input and output files
----------------------

MDCORE takes no arguments and always reads the files `md.dat` and `atoms.dat`.
* `md.dat`:    controls the simulation
* `atoms.dat`: contains the initial configuration, including velocities, charges
               and other quantities that characterize the state of the
               simulation.
`atoms.dat` can be thought of as a restart file.

MDCORE always generates the files `md.out`, `atomsA.out`, `atomsB.out` and
`atoms.out`.
* `md.out`:     reflects the contents of `md.dat`, but amended by the default
                parameters that are not explicitly specified in `md.dat`. This
                is useful to have an overview over all simulation parameters.
                `md.out` also contains annotations that explain each paramter.
                This file is written directly after the code has initialized.
* `atoms*.out`: These files are restart files that have the same format as
                `atoms.dat`. The code writes intermediate alternating restart
                files `atomsA.out` and `atomsB.out`. The frequency of writing
                these files can be controlled by the `file_freq` parameter in
                `md.dat` When the simulation completes, an additional file
                `atoms.out` that contains the final state of the simulation is
                written.

### Structure of `md.dat`

The simulation control file contains two type of elements: key-value pairs
and section. The file always begins with the global section called `Simulation`.

A key-value pair is specified by

    key = "value";

Note that the value has to be enclosed in quotation marks, even if it is a
number and not a string. Lines end with a semicolon (;).

A section starts with the section name followed by `{ };`. Inside the curly
brackets there can be key-value pairs or additional section. As an example:

    OutputNC {
      freq = "10";
    };

Don't forget the semicolon at the end of the section definition.

Comments can be introduced and they start with a `#`.

A full simulation control file (`md.dat`) could look like this:

    Simulation {
      max_time = "100.2";
      OutputNC {
        freq = "10";
      };
    };

### Structure of `atoms.dat`

The configuration file contains sections with different types of data. Each
section starts with the moniker `<---`. The first three sections are fixed,
the following sections are identified by keywords.

Example of a system with four gold atoms:

    <--- Total number of atoms
                    4
    <--- *** The following line is ignored ***
    
    <--- Element, atomic mass, coordinates, group, dissipation, temperature, (next)
      Au     2.04E+04    1.57E+01    1.30E+01    2.60E+00   -2    0    0
      Au     2.04E+04    1.29E+01    1.29E+01    2.42E+00   -2    0    0
      Au     2.04E+04    1.01E+01    1.29E+01    2.43E+00   -2    0    0
      Au     2.04E+04    4.69E-02    1.04E+01    2.43E+00   -2    0    0
      
The columns are explained in the section headers. A `group` is an internal flag
that allows to group atoms together. Groups with negative values are
automatically immovable. 

The `dissipation` and `temperature` columns are parameters for a Langevin
thermostat. These parameters can vary per atom.

The `next` column is only required for molecular systems. It determines the next
atom in this molecule, i.e. can be used to construct a linked-list of atoms
that belong to a certain molecule.

Note that `atomic mass` and `next` are overriden by some MDCORE modules.

Self-documentation
------------------

The code is self-documented. If you misspell a keyword in `md.dat`, then you
will be given a list of all possible keywords with an explanation. If, for
example, you do not know the parameters that the `Langevin` thermostat (see
below) takes, just create a section

    Langevin {
      dummy = "";
    };

and run the code. It will spit out all possible options for the `Langevin`
section.

Simulation flow
---------------

### Code shutdown

The simulation stops at the _time_ specified by `max_time` (in the simulation
control file). Note that this is _not_ the number of time steps, but the time
in the system of units selected. This is particularly useful when using adaptive
time stepping. At this point, the file `atoms.out` and an empty file `DONE`
is written.

MDCORE has a mechanism built-in that allows it to gracefully shut down when the
code hits the wallclock time. In that case, the code _also_ writes the file
`atoms.out`. The existance of the file `atoms.out` does _not_ mean the
simulation has finished. Check for the existence of the `DONE` file to test
whether the simulation has actually hit the time limit.

### Restart

The simulation can be restarted simply by copying `atoms.out` to `atoms.dat`
in a different directory. (All other simulation control files, such as `md.dat`
of course also need to be present in that directory.) Note that the `atoms.out`
file contains the current time. Restarting the simulation then means that the
simulation actually stops at a time `max_time`; `max_time` does _not_ mean the
simulation will run for that duration after startup of the code. This has the
consequence that simulation will not run if it has successfully completed and
is restarted as is.

You can get rid of the time counter by editing the `atoms.out` file and removing
the `time` section. Note that there may be other `time` counters in the file
that should also be removed. For example, output modules (see below) keep track
of the time through their own counters.

### Trajectories

Trajectory information is not automatically dumped to a file. You need to
specify an output module. Presently supported are XYZ, PDB, CFG and NC output.
We strongly suggest to use the NC (NetCDF) output module. It creates a NetCDF
file that follows the [AMBER NetCDF convention](http://ambermd.org/netcdf/nctraj.xhtml).

The statement in the simulation control file that activates output via one of
these modules is

    OutputNC {
      freq = "100";
    };

where `freq` specifies how often a file written in units of real time. For
example, if fs time units are selected, this writes a snapshot every 100 fs.

The are other specialized output modules. For example, it is useful to dump an
overview of energies, stresses, and other themodynamic quantities in the system.
This can be achieved with

    OutputEnergy {
      freq = "100";
      average = "yes";
    };
    
The `average` flag enables computation of averages over the time interval rather
than  dump the values of the state of the simulation when the snapshot is taken.
For simulation in which the cell changes `OutputCell` is useful. Finally,
`OutputTime` produces information on estimates of when the simulation finishes
and how much real time is required to simulate a ns in simulation time.

Note that the file names cannot be changed. Trajectories are typically written
in a file named `traj.nc`. (The extension depends on the file type.) Energies
are written to the file `ener.out`.

## Integrators, thermostats and barostats

The simulation carries out by default a NVE run with a Velocity-Verlet
integrator. You can explicitly specify

    Verlet { };

in the input file, but this corresponds to the default setting. (You should
see this section in the output file `md.out` if no integrator is specified
in the input file.)

The only other integrators available are Langevin integrators. There is
`Langevin`, `Langevin1D`, `LocalLangevin` and `LocalLangevin1D`. The `1D`
integrators thermalize only in one Cartesian direction. The `Local` integrators
apply a different per atom temperature and dissipation, as specified in the
`atoms.dat` file. The standard `Langevin` applies the same parameters (supplied
in the simulation control file) to all atoms.

Other thermostats that are implemented are:
* `BerendsenT`: Berendsen thermostat
* `PetersT`:    Peters thermostat that implements the DPD equations of motion
* `SlidingT`:   A Langevin thermostat that removes the relative velocity of two
                rigidly sliding boundaries. Has to be used in conjuction with
                `SlidingP` (see below).

In addition to thermostats, the following barostats exist:
* `AndersenP`:  Andersen barostat
* `BerendsenP`: Berendsen barostat
* `SlidingP`:   Pressure control for sliding friction simulations

Structure optimization
----------------------

MDCORE presently only implements FIRE. The simulation control file syntax for
structure optimization is

    FIRE {
      fmax = "0.01";
    };

where `fmax` is the convergence criterion that is a threshold on the maximum
residual force. There are other parameters that let you fine tune FIRE's
behavior but the default settings should be appropriate for most purposes.