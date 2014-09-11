#!/bin/bash

#Contornos de las lineas Nodales  para datos de funciones reales en el plano
#Usese seguido del nombre del archivo.
#para datos en Matrix Data. Estoy abusando de que ya se que son los ejes aqui.


Nivel=$1
Name='CentrosPuntos-'$Nivel'.dat'
Energia=$(sed -n $1p EnergiasNelson-0.1875.dat) 
NameFinal='CuerdasPuntos-'$Nivel'.dat'
NamePromedios='PromediosPuntos-'$Nivel'.dat'

echo $Energia
JustCuerdasPara01.x $Name
mv CuerdasParelel.dat $NameFinal
mv PromediosParalel.dat $NamePromedios

echo "Venganza!"


