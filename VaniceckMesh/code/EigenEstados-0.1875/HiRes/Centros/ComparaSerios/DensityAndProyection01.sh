#!/bin/bash

#Contornos de las lineas Nodales  para datos de funciones reales en el plano
#Usese seguido del nombre del archivo.
#para datos en Matrix Data. Estoy abusando de que ya se que son los ejes aqui.

Nivel=$1
Energia=$(sed -n $1p EnergiasNelson-0.1875.dat) 
NameDataDensity=EstadoInterpolado$1-Centros.dat
NameDataProy=EstadoInterpolado$1-PsienQy.dat
Name=`basename $NameDataDensity .dat`
NameTex=$Name"WigneryProy.png"
Using1='(w($1)):(h($2)):3 matrix'
Using2='(w($0)):($2**2)'

gnuplot<<EOF
w(x)=-1.01+10.0*x/1002
h(x)=-1.31+2.62*x/1002
set contour
unset surface
set view map
set cntrparam levels discrete 0.0000
set table "Zeros.dat"
splot "$NameDataDensity" using $Using1
unset table


set xl "q_y"
set yl "p_y"
set palette defined (0 "#000044", 1 "#99FFFF", 2 "white", 3 "red", 4 "#550000")
set size ratio -1
set xr[-1:9]
set yr[-1.3:1.3]
set cbr[-.004:0.004]
set key center  right top

set title "Energia=$Energia, Nivel=$Nivel"
set pointsize 0.4
set term pngcairo enhanced solid lw 1 size 1100,400
set out "$NameTex"
plot "$NameDataDensity" usi $Using1 w image notitle ,  "Zeros.dat"  w l lt 7 lw 0.5 notitle, "$NameDataProy" usi $Using2 w l lc rgb "#219909" lw 2 t "Densidad de Probabilidad"
set out
EOF


echo "Morte!"


