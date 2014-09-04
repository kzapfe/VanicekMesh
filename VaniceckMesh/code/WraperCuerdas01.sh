#!/bin/bash

#Chingadera para ahorrarnos trabajo.

rm Test.dat
NameData=$1
cp $NameData Test.dat
Name=`basename $NameData .dat`
NameSalidaData=$Name"-Cuerdas.dat"

echo $NameSalidaData
echo $NameSalidapng

ParaCuerdasForInterpol01.x Test.dat

mv CuerdasEnMu.dat $NameSalidaData

