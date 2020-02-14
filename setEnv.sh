#!/bin/bash

QTDIR=/home/cfmcginn/Projects/QT/CSAtlasHI
if [[ -d $QTDIR ]]
then
    echo "QTDIR set to '$QTDIR'; if wrong please fix"
    export QTDIR=$QTDIR
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$QTDIR/lib
else
    echo "QTDIR given, '$QTDIR' not found!!! Please fix" 
fi

fastjetStr=$(fastjet-config --cxxflags)
ISFASTJETGOOD=$(echo $fastjetStr | grep "\-I" | wc -l)

fjContribCount=1
if [ -z "$FJCONTRIB__HOME" ];
then
    fjContribCount=0
fi

#Set below in system specific way to your install; if on lxplus and lcgenv fastjet is setup, this would correspond to running:
#setupATLAS
#lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt fastjet" #FASTJET
#And path:
#fastjet-config --prefix
#If not on lxplus, running:
#fastjet-config --prefix
#Will still return path that LOCALFJPATH should be set to below if fastjet is in fact installed and in path

LOCALFJPATH=/home/cfmcginn/Packages/FastJet/fastjet-install
ISLXPORACF=1

if [[ $ISLXPORACF -eq 1 ]]
then
    setupATLAS
    lsetup "root 6.18.04-x86_64-centos7-gcc8-opt" #root6 distribution
    lsetup "root 6.18.04-x86_64-centos7-gcc8-opt" #root6 distribution    
fi

if [[ $ISFASTJETGOOD -eq 0 ]]
then
    if [[ $ISLXPORACF -eq 0 ]]
    then
	if [[ -d $LOCALFJPATH ]]
	then
            export PATH=$PATH:$LOCALFJPATH/bin
	    
	    fastjetStr=$(fastjet-config --cxxflags)
	    ISFASTJETGOOD=$(echo $fastjetStr | grep "\-I" | wc -l)
	else
	    echo "Given '$LOCALFJPATH' is not a valid directory. Please fix"
	fi	
    else
	lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt fastjet"
	fastjetStr=$(fastjet-config --cxxflags)
	ISFASTJETGOOD=$(echo $fastjetStr | grep "\-I" | wc -l)
    fi
    
    if [[ $ISFASTJETGOOD -eq 0 ]]
    then
	echo "FastJet not currently in system path, AND '$LOCALFJPATH' provided does not lead to fastjet-config needed. please fix"
    fi    
fi

if [[ $ISFASTJETGOOD -eq 1 ]]
then
    if [[ $fjContribCount -eq 0 ]]
    then
	echo "FJCONTRIB__HOME is not set"
	echo "ISLXPORACF is set to '$ISLXPORACF' (IF WRONG PLEASE FIX)"
	
	if [[ $ISLXPORACF -eq 0 ]]
	then
	    echo "Setting to export FJCONTRIB__HOME=`fastjet-config --cxxflags`"
	    export FJCONTRIB__HOME=`fastjet-config --cxxflags`
	else
	    lsetup "lcgenv -p LCG_96b x86_64-centos7-gcc8-opt fjcontrib"
	fi	    
    fi
fi

