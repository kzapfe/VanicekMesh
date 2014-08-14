#Version Standalone del notebook llamado "NelsonPotential a la Hernando et Eduardo"
#Probemos usar argumentos y asi:

import sys
import numpy as np
import numpy.linalg as la

datotimestep=sys.argv[1]
#Parametros
N=80 #Puntos por lado en la malla
xmin=-5.
xmax=5.
ymin=-1.
ymax=9.
dt=float(datotimestep)
hbar=0.05
w1=0.1  #el omega de qx^2

#Potencial Nelson
def PotNelson(x,y):
    result=((w1*x*x)/2.+(y-x*x/2.)**2)
    return result


PuntosQ=np.zeros((N*N,2), float)  #Array vacio para los puntos.

#Creamos la malla
cuenta=0
for j in range(0, N):
    for k in range(0, N):
        PuntosQ[cuenta,0]=xmin+(xmax-xmin)*(j+1./2.)/N
        PuntosQ[cuenta,1]=ymin+(ymax-ymin)*(k+1./2.)/N
        cuenta+=1


#Creamos el propagador
Propagator=np.zeros((N*N,N*N), float)
for j in range(0,N*N):
    for k in range(0,N*N):
        Propagator[j,k]=(np.exp(-(la.norm(PuntosQ[j]-PuntosQ[k])**2/(2*dt*hbar))
                                -(PotNelson(PuntosQ[j,0],PuntosQ[j,1])+
                                PotNelson(PuntosQ[k,0],PuntosQ[k,1]))*dt/(hbar*2.0)))

#Normalizamos el propagador
ZD=(xmax-xmin)*(ymax-ymin)
Propagator=Propagator*ZD/(N*N*2.*hbar*dt*np.pi)

#Resolver el Eigensistema
Chunflas=la.eig(Propagator)


#Hay que ver como putos los ordena
explambdaenergias=Chunflas[0].real
energias=np.sort(np.log(explambdaenergias)*hbar/(-dt))
EigenFunciones=np.append(PuntosQ, Chunflas[1].real,1)


NombreEnergia="EnergiasNelson-"+datotimestep+".dat"
NombreEstados="EigenEstadosNelson-"+datotimestep+".dat"
#Guardamos EigenEnergias y EigenVectores
np.savetxt(NombreEnergia, energias)
np.savetxt(NombreEstados, EigenFunciones)
