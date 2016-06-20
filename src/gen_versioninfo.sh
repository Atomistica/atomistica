#! /bin/base

atomistica_revision=$( cd $1 ; git describe --always --tags --dirty )
atomistica_date=$( cd $1 ; git log -1 --format=%ct | awk '{print strftime("%d%b%y",$1)}' )
atomistica_url=$( cd $1 ; git config --get remote.origin.url )
h=`hostname`
m=`uname -m`


if [[ "$3" == "bgxlf_r" ]]; then
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


if [[ $atomistica_revision == "exported" ]]; then
    if [[ -a $1/REV ]]; then
	atomistica_revsion="$atomistica_revision, `cat $1/REV`"
    fi
else
    echo $atomistica_revision > $1/REV
    echo $atomistica_date > $1/REV_DATE
fi

mkdir -p $2
cat<<EOF > $2/versioninfo.f90
module versioninfo
implicit none
integer, private, parameter :: MAXSTRLEN = 1000
character(MAXSTRLEN)  :: atomistica_revision  = "$atomistica_revision"
character(MAXSTRLEN)  :: atomistica_date      = "$atomistica_date"
character(MAXSTRLEN)  :: atomistica_url       = "$atomistica_url"
character(MAXSTRLEN)  :: builddate            = __DATE__ // " " // __TIME__
character(MAXSTRLEN)  :: buildhost            = "$h"
character(MAXSTRLEN)  :: arch                 = "$m"
character(MAXSTRLEN)  :: compileroptions      = "$fortopts"
character(MAXSTRLEN)  :: compilerversion      = "$fortvers"
endmodule versioninfo
EOF
