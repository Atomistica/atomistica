#!/bin/bash
# Rebuild and reinstall atomistica during development using uv
# This script builds a wheel and installs it with force-reinstall

set -e  # Exit on error

echo "Building atomistica wheel with uv..."
uv run --no-sync python -m build --no-isolation -w

echo "Installing atomistica with uv..."
uv pip install dist/atomistica-*.whl --force-reinstall

echo "✓ Successfully rebuilt and installed atomistica"
echo ""
echo "Test with: uv run python -c 'import atomistica; print(\"Success\")'"
