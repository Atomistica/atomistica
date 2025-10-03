#!/bin/bash
# Rebuild and reinstall atomistica during development using uv
# This script builds a wheel and installs it with force-reinstall

set -e  # Exit on error

echo "Ensuring build dependencies are installed..."
uv pip install --quiet build meson-python meson ninja setuptools setuptools-scm

echo "Building atomistica wheel with uv..."
# Use absolute path and run from /tmp to avoid import confusion with build/ directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_PYTHON="$SCRIPT_DIR/.venv/bin/python"
(cd /tmp && "$VENV_PYTHON" -m build --no-isolation -w "$SCRIPT_DIR")

echo "Installing atomistica with uv..."
uv pip install dist/atomistica-*.whl --force-reinstall

echo "✓ Successfully rebuilt and installed atomistica"
echo ""
echo "Test with: uv run python -c 'import atomistica; print(\"Success\")'"
