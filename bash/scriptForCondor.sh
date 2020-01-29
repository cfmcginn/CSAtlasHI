#!/bin/bash

id=0
if [[ $# -eq 1 ]]
then
    id=$1
fi

inputFile=input/rucioJZAllFiles.txt

#source /usatlas/u/cfmcginn/.setEnv.sh

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
export RUCIO_ACCOUNT=cmcginn
lsetup "rucio"
voms-proxy-init --pwstdin <<< 'ZgXAj)C[?I&O.n!M/c3+{_0\;<-l#r9V$GAe'

tar -xzvf inTar.tar.gz

files=()
while read -r line
do
    files+=($line)
done < $inputFile

file=${files[$id]}

if [[ -f $file ]]
then
    echo "RUNNING OVER LOCALLY AVAILABLE FILE $file"
else
    echo "DOWNLOADING FILE $file"
    rucio download $file --no-subdir

    if [[ -f $file ]]
    then
	dummy=0
    else
	echo "DOWNLOAD FAILED, exit 1."
	exit 1
    fi
fi
    
#lsetup "git"
lsetup "root 6.18.04-x86_64-centos7-gcc8-opt"

make clean
make

jzStr=""
if [[ $file == *"19989634"* ]]
then
    jzStr="JZ1"
elif [[ $file == *"19947473"* ]]
then
    jzStr="JZ2"
elif [[ $file == *"19989702"* ]]
then
    jzStr="JZ3"
fi

./bin/clusterToCS.exe $file $file "trk" $jzStr
cp output/*/*.root .

rm $file
