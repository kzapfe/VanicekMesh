//Metodo de 2013 Vanicek et al para hacer super buenas 
//aproximaciones a estados cuanticos.
//Primer intento de hacer el propagador

//primero, todos los include no hechos por mi
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <gsl/gsl_rng.h> //Inicializa el gsl random
#include <gsl/gsl_randist.h> //las distribuciones de gsl
#include <armadillo> //intentemos primero con arma


//#include <ctime> //para ver cuanto se tarda en cada punto
//#include <omp.h> //paralel

// despues todos los include de mi creacion

#include "QuantumConstants.hpp"
#include "BoltzmanforNelson.hpp"
#include "PropagatorFreeParticle.hpp"


using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

  const int PuntosMalla=1000;
 
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);


  cout<<"Hola"<<endl;

  double x,y;
  double NPot;
  double dado;
  double PDFBoltz;
  int contando=0;
  mat propagator=eye<mat>(PuntosMalla,PuntosMalla);
  mat PuntosQ(PuntosMalla,2);
  double factorpropagator=(1./(2.*pi*hbar*dt));
  
  propagator=factorpropagator*propagator;

  

  while(contando<PuntosMalla){
    x=gsl_ran_flat( r, -7, 7);
    y=gsl_ran_flat( r, -2,14);
    dado=gsl_ran_flat( r, 0, 1.0);
    PDFBoltz=Boltzmanpdf(x,y);
    if(PDFBoltz>dado){
      PuntosQ(contando,0)=x;
      PuntosQ(contando,1)=y;
      contando++;
    }
  };

  for(int i=0; i<PuntosMalla-1;i++){
    for(int j=i+1; j<PuntosMalla;j++){  
      propagator(i,j)=PropagatorFreeImagTime(PuntosQ(i,0), PuntosQ(j,0),dt)*
	PropagatorFreeImagTime(PuntosQ(i,1), PuntosQ(j,1),dt);
      propagator(j,i)=propagator(i,j);
    }
  }

  //  cout<<propagator;

  PuntosQ.save("TestPuntos.dat", arma_ascii);

  mat eigenestados;
  vec energias;

  eig_sym(energias, eigenestados, propagator);

  energias.save("LamdaExpEnergias.dat", arma_ascii);

  energias=energias/(PuntosMalla/Z_V);

  energias=log(energias);

  energias=energias/(-dt/hbar);
  
  energias.save("Energias.dat", arma_ascii);
  
  eigenestados.save("EigenStates.dat", raw_ascii);

  return 0;

}

