# Building Atomistica Python Interface

The following instructions explain how to compile, use and test the Python interface to Atomistica.

## Requirements

* Python 3.8 or greater (Python 3.12+ recommended)
* NumPy >= 1.21.0 (NumPy 2.x supported)
* Meson >= 1.1.0
* meson-python >= 0.15.0
* A Fortran compiler (gfortran or ifort)
* A C compiler
* A C++ compiler
* LAPACK library

## Quick Start

### 1. Install build dependencies

Using pip:
```bash
pip install meson-python meson ninja numpy ase
```

Or using uv (recommended):
```bash
uv pip install meson-python meson ninja numpy ase
```

### 2. Build the package

To build a wheel:
```bash
python -m build --no-isolation -w
```

To build and install in development mode:
```bash
pip install -e . --no-build-isolation
```

### 3. Install the built wheel

```bash
pip install dist/atomistica-*.whl
```

### 4. Test the installation

Import Atomistica to verify it works:
```bash
python -c "import atomistica; print('Successfully imported atomistica')"
```

### 5. Run the test suite

```bash
cd tests
python run_tests.py
```

Each test can also be run directly for more diagnostic output:
```bash
python bulk_properties.py
python forces_and_virial.py
python rebo2_molecules.py
```

## Build System Details

Atomistica now uses the Meson build system (via meson-python) instead of the deprecated numpy.distutils. The build configuration is defined in:

* `pyproject.toml` - Python package metadata and build backend configuration
* `meson.build` - Meson build configuration
* `build_helpers/generate_factories.py` - Factory code generation script

## Compiler Configuration

The build system will automatically detect your compilers. If you need to specify particular compilers, you can set environment variables:

For GNU compilers:
```bash
export CC=gcc
export CXX=g++
export FC=gfortran
```

For Intel compilers:
```bash
export CC=icc
export CXX=icpc
export FC=ifort
```

## LAPACK Library

The build system automatically detects LAPACK using pkg-config. If LAPACK is not found automatically, you may need to install it or specify its location:

On Ubuntu/Debian:
```bash
sudo apt-get install liblapack-dev
```

On macOS with Homebrew:
```bash
brew install lapack
```

On Fedora/RHEL:
```bash
sudo dnf install lapack-devel
```

## Troubleshooting

### Build fails with "command not found: meson"

Install Meson and ninja:
```bash
pip install meson ninja
```

### Import error: "symbol not found"

This typically indicates a linking issue. Make sure all dependencies (LAPACK, compiler runtime libraries) are properly installed and accessible.

### NumPy compatibility issues

Atomistica supports both NumPy 1.x and 2.x through a compatibility layer. If you encounter NumPy-related errors, try upgrading to NumPy 2.x:
```bash
pip install 'numpy>=2.0'
```

## Migration from numpy.distutils

If you previously built Atomistica with `setup.py`, note that:

* `setup.cfg` is no longer used for compiler configuration
* The build system now uses `pyproject.toml` and `meson.build`
* You no longer need to run `python setup.py clean` between builds
* OpenMP and other compiler flags are automatically configured by Meson

## Advanced Build Options

To pass additional Fortran flags:
```bash
FFLAGS="-fopenmp" pip install -e . --no-build-isolation
```

To enable verbose build output:
```bash
python -m build --no-isolation -w -v
```

For development, you can use meson directly:
```bash
meson setup builddir
meson compile -C builddir
```
