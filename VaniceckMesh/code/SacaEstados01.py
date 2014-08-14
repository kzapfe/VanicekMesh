# Vamos a graficar, a partir de los datos sacadoas con el programa
# a la Hernando, varios estados chidamente interpolados.
# Besos para Fabricio.

import sys

longitud=len(sys.argv)
#print(longitud)
if (longitud==1):
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

Energias=np.loadtxt("EnergiasNelson0185.dat")
print("Esta es la energia del nivel que pediste")
print(Energias[numerico-1]) #Recuerda que numpy numera desde cero.

Estado=np.loadtxt("EigenEstadosNelson0185.dat", usecols=(0,1,numerico+1))

#Las coordenadas del problema
xx = np.arange(-4.9375, 5, 0.125)
yy =np.arange(-0.9375, 9, 0.125)

#Ponemos la funcion como una matriz.
zz=np.reshape(Estado[:,2], (80,80))
zz=zz.T


from scipy import interpolate

print("Interpolando Datos")
f = interpolate.interp2d(xx, yy, zz, kind='cubic', fill_value=0) 
#las comas, pendejo, no olvidar las putas comas

xnew = np.arange(-5.01, 5.01, 0.01)
ynew = np.arange(-1.01, 9.01, 0.01)
#Nueva malla, mejor interpolda
print("calculando valor de la funcion interpolada")
znew = f(xnew, ynew)

NomFinale="EstadoInterpolado"+numerador+".dat"

np.savetxt(NomFinale,znew)
