# Configuration file for the GNU Compiler Collection (gcc/gfortran).
# Rename to setup.cfg.

[config_fc]
fcompiler=gfortran
f90flags=-cpp -fPIC -ffree-form -ffree-line-length-none -x f95-cpp-input
f77flags=-cpp -fPIC -x f77-cpp-input

[build_ext]
libraries=gfortran

# See the docstring in versioneer.py for instructions. Note that you must
# re-run 'versioneer.py setup' after changing this section, and commit the
# resulting files.

[versioneer]
VCS = git
style = pep440
versionfile_source = src/python/atomistica/_version.py
versionfile_build = atomistica/_version.py
tag_prefix =
