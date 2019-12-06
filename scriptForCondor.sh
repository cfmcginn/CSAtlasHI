#!/bin/bash

id=0
if [[ $# -eq 1 ]]
then
    id=$1
fi

#source /usatlas/u/cfmcginn/.setEnv.sh

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
export RUCIO_ACCOUNT=cmcginn
lsetup "rucio"
lsetup "git"
lsetup "root 6.18.04-x86_64-centos7-gcc8-opt"
voms-proxy-init --pwstdin <<< 'ZgXAj)C[?I&O.n!M/c3+{_0\;<-l#r9V$GAe'

tar -xzvf inTar.tar.gz

make clean
make

files=()

while read -r line
do
    files+=($line)
done < input/rucioJZ2Files.txt

if [[ -f ${files[$id]} ]]
then
    echo "RUNNING OVER LOCALLY AVAILABLE FILE ${files[$id]}"
else
    echo "DOWNLOADING FILE ${files[$id]}"
    rucio download ${files[$id]} --no-subdir
fi
    

./bin/clusterToCS.exe ${files[$id]} ${files[$id]} "trk"
cp output/*/*.root .

rm ${files[$id]}
