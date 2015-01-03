# Fix OS X linking issues - thanks to James Kermode

# recompile, forcing a full rebuild
python setup.py build -f

# see what is included - for me, coulomb_factory_c.o and coulomb_factory_f90.o are missing..
nm build/temp.macosx-10.10-x86_64-2.7/libatomisticalib.a | grep factory

# rebuid static library, including all object files
libtool -static -o build/temp.macosx-10.10-x86_64-2.7/libatomisticalib.a $(find build -name \*.o)

# check results
nm build/temp.macosx-10.10-x86_64-2.7/libatomisticalib.a | grep factory

# finally, relink the extension module (.so)
/usr/bin/clang++ -bundle -undefined dynamic_lookup -L/opt/local/lib -Wl,-headerpad_max_install_names -L/opt/local/lib/db48 build/temp.macosx-10.10-x86_64-2.7/Users/jameskermode/Code/atomistica/src/python/c/py_f.o build/temp.macosx-10.10-x86_64-2.7/Users/jameskermode/Code/atomistica/src/python/c/particles.o build/temp.macosx-10.10-x86_64-2.7/Users/jameskermode/Code/atomistica/src/python/c/neighbors.o build/temp.macosx-10.10-x86_64-2.7/Users/jameskermode/Code/atomistica/src/python/c/coulomb.o build/temp.macosx-10.10-x86_64-2.7/Users/jameskermode/Code/atomistica/src/python/c/coulomb_callback.o build/temp.macosx-10.10-x86_64-2.7/Users/jameskermode/Code/atomistica/src/python/c/potential.o build/temp.macosx-10.10-x86_64-2.7/Users/jameskermode/Code/atomistica/src/python/c/analysis.o build/temp.macosx-10.10-x86_64-2.7/Users/jameskermode/Code/atomistica/src/python/c/atomisticamodule.o -Lbuild/temp.macosx-10.10-x86_64-2.7 -latomisticalib -lblas -llapack -o build/lib.macosx-10.10-x86_64-2.7/_atomistica.so -L/opt/local/lib/gcc48 -lgfortran
