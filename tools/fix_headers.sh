#! /bin/sh

ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

DELIM="======================================================================"

# .c/.h files
for i in `find src/lammps src/potentials src/python src/support -name "*.c"` `find src/lammps src/potentials src/python src/support -name "*.h"`; do

    if [ "$(grep $DELIM $i | wc -l)" -gt 0 ]; then

    sed "/\/\* ${DELIM}/,/${DELIM} \*\//d" $i > $i.tmp
    cat ${ROOT}/c_header.txt $i.tmp > $i
    rm $i.tmp

    else

    cat ${ROOT}/c_header.txt $i > $i.tmp
    mv $i.tmp $i

    fi

done

# .f90 files
for i in `find src/lammps src/potentials src/python src/support -name "*.f90"`; do

    if [ "$(grep $DELIM $i | wc -l)" -gt 0 ]; then

    sed "/!! ${DELIM}/,/!! ${DELIM}/d" $i > $i.tmp
    cat ${ROOT}/f_header.txt $i.tmp > $i
    rm $i.tmp

    else

    cat ${ROOT}/f_header.txt $i > $i.tmp
    mv $i.tmp $i

    fi

done
