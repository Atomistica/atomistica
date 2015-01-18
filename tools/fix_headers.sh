#! /bin/sh

ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

DELIM="======================================================================"
NDELIM1="N====================================================================="

# .c/.h files
for i in `find src/core src/lammps src/potentials src/notb src/python src/special src/standalone src/support -name "*.c*"` `find src/lammps src/potentials src/notb src/python src/special src/support -name "*.h"`; do

    if [ "$(grep $DELIM $i | wc -l)" -gt 0 ]; then

    echo "GPL: $i"

    sed "/\/\* ${DELIM}/,/${DELIM} \*\//d" $i > $i.tmp
    cat ${ROOT}/c_header.txt $i.tmp > $i
    rm $i.tmp

    else

    echo "GPL: $i"

    cat ${ROOT}/c_header.txt $i > $i.tmp
    mv $i.tmp $i

    fi

done

# .f90 files
for i in `find src/core src/lammps src/potentials src/notb src/python src/special src/standalone src/support -name "*.f90"`; do

    if [ "$(grep ${NDELIM1} $i | wc -l)" -gt 0 ]; then

    echo "proprietary: $i"

    sed "/!! ${NDELIM1}/,/!! ${DELIM}/d" $i > $i.tmp
    cat ${ROOT}/f_header_n.txt $i.tmp > $i
    rm $i.tmp

    else

    if [ "$(grep ${DELIM} $i | wc -l)" -gt 0 ]; then

    echo "GPL: $i"

    sed "/!! ${DELIM}/,/!! ${DELIM}/d" $i > $i.tmp
    cat ${ROOT}/f_header.txt $i.tmp > $i
    rm $i.tmp

    else

    echo "GPL: $i"

    cat ${ROOT}/f_header.txt $i > $i.tmp
    mv $i.tmp $i

    fi

    fi

done

# .py files
for i in `find src/core src/lammps src/potentials src/notb src/python src/special src/standalone src/support tests -name "*.py"`; do

    if [ "$(grep ${NDELIM1} $i | wc -l)" -gt 0 ]; then

    echo "proprietary: $i"

    sed "/# ${NDELIM1}/,/# ${DELIM}/d" $i > $i.tmp
    cat ${ROOT}/py_header_n.txt $i.tmp > $i
    rm $i.tmp

    else

    if [ "$(grep ${DELIM} $i | wc -l)" -gt 0 ]; then

    echo "GPL: $i"

    sed "/# ${DELIM}/,/# ${DELIM}/d" $i > $i.tmp
    cat ${ROOT}/py_header.txt $i.tmp > $i
    rm $i.tmp

    else

    echo "GPL: $i"

    cat ${ROOT}/py_header.txt $i > $i.tmp
    mv $i.tmp $i

    fi

    fi

done
