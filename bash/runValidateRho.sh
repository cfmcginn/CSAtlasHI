#!/bin/bash

rhoFile=/home/cfmcginn/Samples/QT/rhoFile_HICS_HIEventShapeWeighted_iter0_20200205.root
inFile=/home/cfmcginn/Samples/QT/run.root

./bin/validateRho.exe $rhoFile $inFile
