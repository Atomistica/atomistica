# Development Guide

## Quick Rebuild During Development

When actively developing Atomistica, use these scripts to quickly rebuild and reinstall:

### Using standard pip/venv:
```bash
./rebuild.sh
```

### Using uv:
```bash
./rebuild-uv.sh
```

Both scripts will:
1. Build a new wheel with `python -m build --no-isolation -w`
2. Force-reinstall it with `pip install --force-reinstall`
3. Show a test command to verify the installation

## Why Not Use Editable Installs?

Meson's editable install mode (`pip install -e .`) can have caching issues, especially with the factory code generation and Fortran module dependencies. The wheel rebuild workflow is more reliable for development.

## Full Development Setup

```bash
# Clone the repository
git clone https://github.com/Atomistica/atomistica.git
cd atomistica

# Create a virtual environment (using venv or uv)
python -m venv .venv
source .venv/bin/activate
# OR
uv venv
source .venv/bin/activate

# Install build dependencies
pip install meson-python meson ninja numpy ase build
# OR
uv pip install meson-python meson ninja numpy ase build

# Build and install
./rebuild.sh
# OR
./rebuild-uv.sh

# Run tests
cd tests
python run_tests.py
```

## Making Changes

1. Edit source files (Fortran, C, C++, or Python)
2. Run `./rebuild.sh` or `./rebuild-uv.sh`
3. Test your changes
4. Repeat

## Common Development Tasks

### Running a single test:
```bash
cd tests
python bulk_properties.py
```

### Verbose build output:
```bash
python -m build --no-isolation -w -v
```

### Clean build artifacts:
```bash
rm -rf .mesonpy-* build dist *.egg-info
```

### Using Meson directly:
```bash
meson setup builddir
meson compile -C builddir
```

## Build System Files

- `pyproject.toml` - Python package metadata and build backend configuration
- `meson.build` - Meson build configuration (Fortran/C/C++ compilation)
- `build_helpers/generate_factories.py` - Factory code generation script
- `src/gen_versioninfo.sh` - Version info generation script
