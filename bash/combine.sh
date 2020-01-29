#!/bin/bash

if [[ $# -ne 1 ]]
then
    echo "Usage: bash combine.sh <inDir>. exit 1"
    exit 1
elif [[ -d $1 ]]
then
    count=$(ls $1/spectraComp_*png | wc -l)
    if [[ $count -lt 1 ]]
    then
	echo "Given directory \'$1\' contains no valid png of form \'spectraComp_*png\'. exit 1"
	exit 1
    fi
else
    echo "Given input \'$1\' is not a valid directory. exit 1"
    exit 1
fi

centStr=()
globalStr=()

for i in $1/spectraComp_*png
do
    tempGlobe=${i#*spectraComp_}
    while [[ $tempGlobe == *"_"* ]]
    do
	tempGlobe=${tempGlobe%_*}
    done
	
    tempCent=${i#*_Cent}
    while [[ $tempCent == *"_"* ]]
    do
	tempCent=${tempCent%_*}
    done
    

    globeFound=0
    for i in ${globalStr[@]}
    do
	if [[ $i == $tempGlobe ]]
	then
	    globeFound=1
	    break
	fi
    done

    if [[ $globeFound -eq 0 ]]
    then
	globalStr+=($tempGlobe)
    fi

    centFound=0
    for i in ${centStr[@]}
    do
	if [[ $i == $tempCent ]]
	then
	    centFound=1
	    break
	fi
    done

    if [[ $centFound -eq 0 ]]
    then
	centStr+=($tempCent)
    fi
done


echo "Global: "
for i in ${globalStr[@]}
do
    echo " $i"
done

pos=-1
for i in ${centStr[@]}
do
    pos=$((pos + 1))
done
    
centStr2=()
while [[ $pos -ge 0 ]]
do
    centStr2+=(${centStr[$pos]})
    pos=$((pos - 1))
done

echo "Cent: "
pos=0
for i in ${centStr2[@]}
do
    echo " $i"
    pos=$((pos + 1))
done

for i in ${globalStr[@]}
do
    comboStr=""
    for j in ${centStr2[@]}
    do
	count=$(ls $1/spectraComp_*$i*$j*.png | wc -l)
	file=$(ls $1/spectraComp_*$i*$j*.png)
	if [[ $count -eq 1 ]]
	then
	    comboStr="$comboStr $file"
	else
	    echo "UHOH '$1/*$i*$j' returns not just 1 file. exit 1"
	    exit 1
	fi	
    done

    if [[ $pos -le 4 ]]
    then
	montage $comboStr -tile "$pos"x1 -geometry +0-0 $1/spectraCompCombo_"$i"_"$j".png
    elif [[ $pos -le 6 ]]
    then
	montage $comboStr -tile 3x2 -geometry +0-0 $1/spectraCompCombo_"$i"_"$j".png
    elif [[ $pos -le 8 ]]
    then
	montage $comboStr -tile 4x2 -geometry +0-0 $1/spectraCompCombo_"$i"_"$j".png
    else
	echo "MORE THAN 8 no COMBO AVAIL. exit 1"
	exit 1
    fi
done

rm -f $1/spectraComp_*.png
