#!/bin/bash
# Rebuild and reinstall atomistica during development
# This script builds a wheel and installs it with force-reinstall

set -e  # Exit on error

echo "Ensuring build dependencies are installed..."
pip install --quiet build meson-python meson ninja setuptools setuptools-scm

echo "Building atomistica wheel..."
# Use absolute path and run from /tmp to avoid import confusion with build/ directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
(cd /tmp && python -m build --no-isolation -w "$SCRIPT_DIR")

echo "Installing atomistica..."
pip install dist/atomistica-*.whl --force-reinstall

echo "✓ Successfully rebuilt and installed atomistica (Fortran) + atomistica_cpp (C++)"
echo ""
echo "Test Fortran:   python -c 'import atomistica; print(\"OK\")'"
echo "Test C++:       python -c 'import atomistica_cpp; print(\"OK\")'"
echo "Or with venv:   .venv/bin/python -c 'import atomistica; import atomistica_cpp; print(\"OK\")'"
