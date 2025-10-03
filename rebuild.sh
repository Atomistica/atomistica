#!/bin/bash
# Rebuild and reinstall atomistica during development
# This script builds a wheel and installs it with force-reinstall

set -e  # Exit on error

echo "Building atomistica wheel..."
python -m build --no-isolation -w

echo "Installing atomistica..."
pip install dist/atomistica-*.whl --force-reinstall

echo "✓ Successfully rebuilt and installed atomistica"
echo ""
echo "Test with: python -c 'import atomistica; print(\"Success\")'"
