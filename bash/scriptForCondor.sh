#!/bin/bash

id=0
if [[ $# -eq 1 ]]
then
    id=$1
fi

inputFile=input/rucioJZAllFiles_TEST_20200604.txt

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

echo "RUN LS A"
ls
echo "RUN LS A1"
ls bin/
#./bin/clusterToCS.exe $file2 $file2 "trk" $jzStr
sed -i -e "s@REPINFILENAME@$file2@g" input/defaultConfig.txt
sed -i -e "s@PROCESS@$1@g" input/defaultConfig.txt
./bin/makeClusterTree.exe input/defaultConfig.txt
echo "RUN LS B"
ls
cp output/*/*.root .


echo "RUN LS C"
rm $file2

