#!/bin/bash

id=0
if [[ $# -eq 1 ]]
then
    id=$1
fi

inputFile=input/rucioJZAllFiles_20200604.txt

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
file2=${file#*:}

if [[ -f $file ]]
then
    echo "RUNNING OVER LOCALLY AVAILABLE FILE $file"
    file2=$file
else
    echo "DOWNLOADING FILE $file"
    rucio download $file --no-subdir

    echo " Checking \'$file2\'"
    if [[ -f $file2 ]]
    then       
	dummy=0
    else
	echo "DOWNLOAD FAILED, exit 1."
	rm -f *.root
	exit 1
    fi
fi
    
#lsetup "git"
lsetup "root 6.18.04-x86_64-centos7-gcc8-opt"
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt fastjet"
lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt fjcontrib"
export QTDIR=$PWD
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD:$PWD/lib

make clean
make

#./bin/clusterToCS.exe $file2 $file2 "trk" $jzStr
sed -i -e "s@REPINFILENAME@$file2@g" input/makeClusterTree.config
sed -i -e "s@PROCESS@$1@g" input/makeClusterTree.config
./bin/makeClusterTree.exe input/makeClusterTree.config
echo "RUN LS B"
ls
cp output/*/*.root .


echo "RUN LS C"
rm $file2

