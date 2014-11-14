"""Es mas facil si parseamos los datos a una matriz cuadrada para armadillo.
Este script agarra la lista de datos que hizo Eduardo a la gnuplot, estilo
ya sabes, (x,y,f(x,y))  y la convierte en la matriz cuadrada con
ejes x et  y"""

import sys

longitud=len(sys.argv)
if (longitud<2):
    print('Tienes que dar un nombre de archivo de datos')
    sys.exit()

nombre=sys.argv[1]

print(nombre)

import numpy as np


#Las coordenadas del problema
xx = np.loadtxt(nombre, usecols=([0]))
yy =np.loadtxt(nombre, usecols=([1]))
Estado=np.loadtxt(nombre, usecols=([2]))
##Numpy ignora automaticamente los espacios en blanco, asi que pues ya nos
##ahorramos un script de limpiar datos.

NumP=int(np.sqrt(xx.shape[0]))

Estado=np.reshape(Estado, (NumP,NumP))
Estado=Estado.T

outname=nombre.split(".")[0]
outname=outname+"-valores.dat"

np.savetxt(outname, Estado)

#Esto bastaba hacerlo una sola vez
#xlimpio=np.unique(xx)
#ylimpio=np.unique(yy)

#coordenadas=np.empty([2,NumP])
#coordenadas[0]=xlimpio
#coordenadas[1]=ylimpio

#coordenadas=np.transpose(coordenadas)

#np.savetxt("Coordenadas01.dat", coordenadas)


#print(xlimpio)
