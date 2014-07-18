//la distribucion de Boltzman para el potencial de Nelson
#ifndef __PropagatorFreeParticle__
#define __PropagatorFreeParticle__

//necesita QuantumConstants

using namespace std;
using namespace arma;

double PropagatorFreeImagTime(double x, double xprima, double deltatau){
  //Este es el propagador con tiempo Imaginario. 
  //En realidad es un numero real, que depende de puros numeros reales
  //acorde al ITNUM de Vanicek y Hernando
  //masa=1 de momento
  //esta parte es sin el factor Z_V/N que incluye el articulo
  double prefactor;
  double result;
  prefactor=1.0/(2.*pi*hbar*deltatau);
  result=sqrt(prefactor)*exp(-prefactor*pi*(x-xprima)*(x-xprima));
  return result;
  
  
}

double UNormPropagatorFreeImagTime(rowvec q, rowvec qprima, double deltatau){
  //Este es el propagador con tiempo Imaginario. 
  //En realidad es un numero real, que depende de puros numeros reales
  //acorde al ITNUM de Vanicek y Hernando
  //masa=1 de momento
  //solo el exponente, hbar=1, m=1
  double result;
  result=exp(-(norm(q-qprima)*norm(q-qprima))/(2.*deltatau*hbar));
  return result;  
  
}


mat FreePropagatorImagTime(mat puntos, double deltatau){
  //Maybe this wont work.
  int tamanho=puntos.n_rows;
  mat result=zeros<mat>(tamanho,tamanho);
 
  return result;
  
}

double UNormPropOAITime(double x, double xprima, 
			double y, double yprima, 
			double deltatau){
  
  double result;
  result=exp((-PotArm(x,y)-PotArm(xprima,yprima))*deltatau/(2.0*hbar));
  return result;  
  
}

double PropOA1D(double x, double xprima, 
	    	double deltatau){
  
  double w=1.0;
  double result;
  result=exp(-(x*x+xprima*xprima)*(w/2.)*deltatau/(2.0*hbar));
  return result;  
  
}




#endif
