#!/bin/bash

fileCalo=/home/cfmcginn/MITHIG/tempPlots2019/Dec2plots/mc16_5TeV.root
fileATLAS=/home/cfmcginn/MITHIG/tempPlots2019/Dec2plots/mc16_5TeV.root

./bin/clusterToCS.exe $fileCalo $fileATLAS "trk"
