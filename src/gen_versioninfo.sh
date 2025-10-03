#! /bin/bash

# Get version from setuptools_scm
atomistica_revision=$( cd $1/.. ; python3 -c "try:
    from setuptools_scm import get_version
    print(get_version(root='..', relative_to=__file__))
except:
    print('0.0.0')
" 2>/dev/null || echo "0.0.0" )

# Get current date for atomistica_date
atomistica_date=$(date '+%Y-%m-%d')
atomistica_url=$( cd $1 ; git config --get remote.origin.url )
if [ -z "$atomistica_url" ]; then
  atomistica_url="N/A"
fi
h=`hostname`
m=`uname -m`


if [ "$3" = "bgxlf_r" ]; then
    fortvers1=`$3 -qversion | head -n 1 | tail -n 1`
    fortvers2=`$3 -qversion | head -n 2 | tail -n 1`
    fortvers="$fortvers1; $fortvers2"
elif [ -z "$3" ] || [ "$3" = "Python" ]; then
    fortvers="meson-python"
else
    fortvers=`$3 --version | head -n 1 2>/dev/null || echo "unknown"`
fi


fortopts=""
n=0
for i in "$@"; do
   let n=$n+1

   if [ $n -gt 3 ]; then
      fortopts="$fortopts$i "
   fi
done


mkdir -p $2

# Get current date and time
builddate=$(date '+%b %d %Y %H:%M:%S')

cat<<EOF > $2/versioninfo.f90
module versioninfo
implicit none
integer, private, parameter :: MAXSTRLEN = 1000
character(MAXSTRLEN)  :: atomistica_revision  = "$atomistica_revision"
character(MAXSTRLEN)  :: atomistica_date      = "$atomistica_date"
character(MAXSTRLEN)  :: atomistica_url       = "$atomistica_url"
character(MAXSTRLEN)  :: builddate            = "$builddate"
character(MAXSTRLEN)  :: buildhost            = "$h"
character(MAXSTRLEN)  :: arch                 = "$m"
character(MAXSTRLEN)  :: compileroptions      = "$fortopts"
character(MAXSTRLEN)  :: compilerversion      = "$fortvers"
endmodule versioninfo
EOF
