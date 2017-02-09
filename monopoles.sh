#! /bin/bash
cluster=$1
process=$2
basedir=$3

mkdir $basedir/beta.05
cd $basedir

monopoles 1000 5e-5 .05 beta.05/monopoles_$process.root >& /dev/null

exit 0
