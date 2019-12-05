#!/bin/bash

id=0
if [[ $# -eq 1 ]]
then
    id=$1
fi

source /usatlas/u/cfmcginn/.setEnv.sh

#setupATLAS
#lsetup "root 6.18.04-x86_64-centos7-gcc8-opt"

tar -xzvf inTar.tar.gz

make clean
make

files=()

while read -r line
do
    files+=($line)
done < input/jz2Files.txt

./bin/clusterToCS.exe ${files[$id]} ${files[$id]} "trk"
cp output/*/*.root .
