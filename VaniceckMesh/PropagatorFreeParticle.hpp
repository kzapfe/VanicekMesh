//la distribucion de Boltzman para el potencial de Nelson
#ifndef __PropagatorFreeParticle__
#define __PropagatorFreeParticle__

//necesita QuantumConstants

using namespace std;

double PropagatorFreeImagTime(double x, double xprima, double deltatau){
  //Este es el propagador con tiempo Imaginario. 
  //En realidad es un numero real, que depende de puros numeros reales
  //acorde al ITNUM de Vanicek y Hernando
  //masa=1 de momento
  //esta parte es sin el factor Z_V/N que incluye el articulo
  double prefactor;
  double result;
  prefactor=1.0/(2.*pi*hbar*dt);
  result=sqrt(prefactor)*exp(-prefactor*pi*(x-xprima)*(x-xprima));
  return result;
  
  
}



#endif
