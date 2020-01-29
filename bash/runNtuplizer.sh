#!/bin/bash

export DOGLOBALDEBUGROOT=0 #SETS AN ENV VAR TO PICKUP GLOBALLY
if [[ $# -eq 1 ]]
then
    export DOGLOBALDEBUGROOT=$1
fi

configFile=/home/cfmcginn/Projects/QT/input/defaultConfig.txt
./bin/makeClusterTree.exe $configFile
