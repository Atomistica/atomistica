#! /bin/bash

mdcore_rev=$( cd $1 ; git describe --always --tags --dirty )
mdcore_url=$( cd $1 ; git config --get remote.origin.url )
libAtoms_rev=`cat $1/libAtoms/REV | grep Revision | awk '{ print $2 }'`
libAtoms_url=`cat $1/libAtoms/REV | grep URL | awk '{ print $2 }'`
h=`hostname`
m=`uname -m`


if [[ "$3" == "xlf_r" ]]; then
    fortvers1=`$3 -qversion | head -n 1 | tail -n 1`
    fortvers2=`$3 -qversion | head -n 2 | tail -n 1`
    fortvers="$fortvers1; $fortvers2"
else if [[ "$3" == "Python" ]]; then
    fortvers="numpy-distutils"
else
    fortvers=`$3 --version | head -n 1`
fi
fi


fortopts=""
n=0
for i in "$@"; do
   let n=$n+1

   if [ $n -gt 3 ]; then
      fortopts="$fortopts$i "
   fi
done


if [[ $mdcore_rev == "exported" ]]; then
    if [[ -a $1/REV ]]; then
	mdcore_rev="$mdcore_rev, `cat $1/REV`"
    fi
else
    echo $mdcore_rev > $1/REV
fi

cat<<EOF > $2/versioninfo.f90
module versioninfo
implicit none
integer, private, parameter :: MAXSTRLEN = 1000
character(MAXSTRLEN)  :: mdcore_revision      = "$mdcore_rev"
character(MAXSTRLEN)  :: mdcore_url           = "$mdcore_url"
character(MAXSTRLEN)  :: libAtoms_revision    = "$libAtoms_rev"
character(MAXSTRLEN)  :: libAtoms_url         = "$libAtoms_url"
character(MAXSTRLEN)  :: builddate            = __DATE__ // " " // __TIME__
character(MAXSTRLEN)  :: buildhost            = "$h"
character(MAXSTRLEN)  :: arch                 = "$m"
character(MAXSTRLEN)  :: compileroptions      = "$fortopts"
character(MAXSTRLEN)  :: compilerversion      = "$fortvers"
endmodule versioninfo
EOF
