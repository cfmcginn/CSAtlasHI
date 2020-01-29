#!/bin/bash

fileCalo=/atlasgpfs01/usatlas/data/cfmcginn/CaloCluster/user.cmcginn.jetreader.20191203.193833.mc16_5TeV.420012.HIJZ2_myOutput.root/user.cmcginn.19947473._000001.myOutput.root
fileATLAS=/atlasgpfs01/usatlas/data/cfmcginn/CaloCluster/user.cmcginn.jetreader.20191203.193833.mc16_5TeV.420012.HIJZ2_myOutput.root/user.cmcginn.19947473._000001.myOutput.root

pos=0
for i in /atlasgpfs01/usatlas/data/cfmcginn/CaloCluster/user.cmcginn.jetreader.20191203.193833.mc16_5TeV.420012.HIJZ2_myOutput.root/user.cmcginn.19947473._00000*.root
do
    ./bin/clusterToCS.exe $i $i "trk" >& trk$pos.log
    pos=$((pos+1))
    exit 1
done
