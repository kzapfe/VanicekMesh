#!/bin/bash

#Contornos de las lineas Nodales  para datos de funciones reales en el plano
#Usese seguido del nombre del archivo.
#para datos en Matrix Data. Estoy abusando de que ya se que son los ejes aqui.


Energia=$1
NameData=EstadoInterpolado$1.dat
Name=`basename $NameData .dat`
NameTex=$Name".png"
Using='(f($1)):(h($2)):3 matrix'
Energia=$(sed -n $1p EnergiasNelson-0.1875.dat) 

echo $Energia

gnuplot<<EOF
#Funciones auxiliares para las coordendas
f(x)=-5.+10.*x/1002
h(x)=-1.+10.*x/1002
omega=0.1
E=$Energia
print E
#Funciones auxiliares para la caustica
ymax(x)=0.5*(x**2+sqrt(4.*E-2.*omega*x**2))
ymin(x)=0.5*(x**2-sqrt(4.*E-2.*omega*x**2))
set samples 500

set xl "q_x"
set yl "q_y"
set palette defined (0 "#000044", 1 "#99FFFF", 2 "white", 3 "red", 4 "#550000")
set size ratio -1
set xr[-5:5]
set yr[-1:9]
set cbr[-.15:0.15]
set grid front

set title "Energia=$Energia, Estado $1"

set pointsize 0.4
set term pngcairo enhanced solid lw 1 size 750,750 font "/usr/share/fonts/TTF/luximr.ttf" 
set out "$NameTex"
plot "$NameData" usi $Using w image notitle, ymax(x) ls 7 lw 2 notitle, ymin(x) ls 7 lw 2 t "Classical Caustic"
set out
EOF


echo "Morte!"


