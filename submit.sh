#!/bin/bash

cFile=input/condorTEMPLATE.job

DATE=`date +%Y%m%d`
TIME=`date +%H%M%S`
cDir=condorDir/condor_"$DATE"_"$TIME"

mkdir -p $cDir

cFile2=$cDir/condor_"$DATE"_"$TIME".job
cp $cFile $cFile2

sed -i -e "s@INITDIR@$PWD/$cDir@g" $cFile2
sed -i -e "s@Queue@@g" $cFile2

files=()

for i in include/*.h
do
    files+=($i)
done

for i in src/*.C
do
    files+=($i)
done

for i in input/*
do
    files+=($i)
done

files+=(Makefile)

inputStr=""
for i in ${files[@]}
do 
    echo $i
    inputStr="$inputStr $i"
done

tar -czvf inTar.tar.gz $inputStr
mv inTar.tar.gz $cDir

count=$(cat input/jz2Files.txt | wc -l)

echo "transfer_input_files = inTar.tar.gz,scriptForCondor.sh" >> $cFile2
echo "Queue $count" >> $cFile2

cp scriptForCondor.sh $cDir
