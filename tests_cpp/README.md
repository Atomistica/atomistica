# atomistica_cpp Python Tests

Python-level tests for the `atomistica_cpp` package.  Tests run against
the installed wheel, so build and install first:

```bash
cd atomistica_cpp
./rebuild.sh            # or ./rebuild-uv.sh
```

## Running

```bash
cd tests_cpp
pip install pytest numpy ase
pytest                  # all tests
pytest test_cpp_forces_and_virial.py   # single file
pytest -k Tersoff       # filter by name
```

## Requirements

- `atomistica_cpp` installed (from wheel)
- `ase >= 3.15`
- `numpy >= 1.21`
- Test data files from the Fortran test suite
  (`../atomistica_fortran/tests/`) for EAM and amorphous carbon tests

## File Structure

| File | Tests |
|---|---|
| `test_cpp_forces_and_virial.py` | Numerical force/stress consistency for all potentials |
| `test_cpp_bulk_properties.py` | Elastic constants and lattice parameters |
| `test_cpp_neighbor_list.py` | NeighborList correctness, PBC, Verlet shell |
| `test_cpp_eam.py` | EAM loading, forces, stress, crash cases |
| `test_cpp_coulomb.py` | Coulomb energy/forces/stress |
| `test_cpp_pbc.py` | PBC force sums, wrapping invariance |

## Notes

- Tests are skipped automatically when ASE is not installed.
- Test data files (`.eam`, `.cfg`) are looked up relative to
  `../atomistica_fortran/tests/`.  Individual tests using these files
  are skipped with a clear message when the data is not found.
- REBO2 forces are intentionally simplified in the C++ implementation
  (angular derivatives omitted).  Force tests for REBO2 use only dimers
  where the simplified forces are exact.
