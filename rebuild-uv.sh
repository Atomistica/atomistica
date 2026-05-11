#!/bin/bash
# Rebuild and reinstall atomistica during development using uv

set -e

echo "Ensuring build dependencies are installed..."
uv pip install --quiet build meson-python meson ninja setuptools setuptools-scm

echo "Building atomistica wheel with uv..."
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_PYTHON="$SCRIPT_DIR/.venv/bin/python"
(cd /tmp && "$VENV_PYTHON" -m build --no-isolation -w "$SCRIPT_DIR")

echo "Installing atomistica with uv..."
uv pip install dist/atomistica-*.whl --force-reinstall

echo "✓ Successfully rebuilt and installed atomistica"
echo ""
echo "Test C++:    .venv/bin/python -c 'import atomistica_cpp; print(\"OK\")'"
echo "Test compat: .venv/bin/python -c 'import atomistica; print(\"OK\")'"
echo ""
echo "Note: Don't use 'uv run' — it will reinstall as editable and fail."
echo "      Always use .venv/bin/python directly."
