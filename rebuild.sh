#!/bin/bash
# Rebuild and reinstall atomistica during development

set -e

echo "Ensuring build dependencies are installed..."
pip install --quiet build meson-python meson ninja setuptools setuptools-scm

echo "Building atomistica wheel..."
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
(cd /tmp && python -m build --no-isolation -w "$SCRIPT_DIR")

echo "Installing atomistica..."
pip install dist/atomistica-*.whl --force-reinstall

echo "✓ Successfully rebuilt and installed atomistica"
echo ""
echo "Test C++:   python -c 'import atomistica_cpp; print(\"OK\")'"
echo "Test compat: python -c 'import atomistica; print(\"OK\")'"
