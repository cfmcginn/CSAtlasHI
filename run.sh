#!/bin/bash

fileCalo=/home/cfmcginn/MITHIG/tempPlots2019/Nov25plots/user.cmcginn.19848289.EXT0.MERGED.caloTree_HIClusterSubtraction.root
fileATLAS=/home/cfmcginn/MITHIG/tempPlots2019/Nov25plots/userATLAS_MERGED.root

./bin/clusterToCS.exe $fileCalo $fileATLAS
