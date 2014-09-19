# Vamos a graficar, a partir de los datos sacadoas con el programa
# a la Hernando, varios estados chidamente interpolados.
# Besos para Fabricio.

import sys

longitud=len(sys.argv)
#print(longitud)
if (longitud<2):
    print('Tienes que dar un argumento numerico entre 1 y 350')
    sys.exit()
    

numerador=sys.argv[1]
numerico=int(numerador)

hbar=0.05
xmin=-5.
xmax=5.
ymin=-1.
ymax=9.

import numpy as np

Energias=np.loadtxt("Energias-0.03-256.dat")
print("Esta es la energia del nivel que pediste")
print(Energias[numerico-1]) #Recuerda que numpy numera desde cero.

Estado=np.loadtxt("Estados-0.03-256.dat", usecols=(0,1,numerico+1))

#Las coordenadas del problema
xx = np.loadtxt("Estados-0.03-256.dat", usecols=([0]))
yy =np.loadtxt("Estados-0.03-256.dat", usecols=([1]))
#Mira, detallito pendejo de python: un array de UN SOLO PUTO numero tiene que ser
#especificao como array-
NumP=int(np.sqrt(xx.shape[0]))

print("Longitud por lado de los datos originales:", NumP)
print(xx[0], xx[NumP*NumP-1])
print(yy[0], yy[NumP*NumP-1])
#Ponemos la funcion como una matriz.
zz=np.reshape(Estado[:,2], (NumP,NumP))
zz=zz.T


from scipy import interpolate

print("Interpolando Datos")
#f = interpolate.interp2d(xx, yy, zz, kind='cubic', fill_value=0)
f2 = interpolate.bisplrep(xx, yy, zz,  s=0.001, eps=0.00001)

#las comas, pendejo, no olvidar las putas comas
print("Terminamos la Interpolacion")

xnew = np.arange(-5.01, 5.01, 0.005)
ynew = np.arange(-1.01, 9.01, 0.005)
#Nueva malla, mejor interpolda0
print("calculando valor de la funcion interpolada")
#znew=f(xnew, ynew)
znew = interpolate.bisplev(xnew, ynew, f2)

NomFinale="EstadoInterpoladoMedioDenso"+numerador+".dat"

np.savetxt(NomFinale,znew)
