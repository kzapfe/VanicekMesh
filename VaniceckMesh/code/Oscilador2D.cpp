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
#include <gsl/gsl_cdf.h> //las distribuciones de gsl
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

  const int Puntos1D=120;
  const int PuntosMalla=Puntos1D*(Puntos1D+1)/2;
  //En este test, vamos a usar hbar=1

  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);


  cout<<"Probando exactamente lo de V et H."<<endl;

  double x, dumy;
  double PDFBoltz;
  int contando=0;
  
  mat propagator=eye<mat>(PuntosMalla,PuntosMalla);
  
  vec PuntosX(Puntos1D), PuntosY(Puntos1D);
  mat PuntosQ(PuntosMalla,2);

  double factorpropagator=sqrt(1./(2.*pi*hbar*dt));
  
  const double  Z_O1D=sqrt(2.*pi/dt);

  for(int i=0; i<Puntos1D;i++){
    //exactamente igual que en los adendums del paper
    //muestreo "perfecto"
    dumy=(double)(i+0.5)/(double)Puntos1D;
    PuntosX(i)=sqrt(1./dt)*
      gsl_cdf_ugaussian_Pinv(dumy);
    cout<<"Llevo "<<i+1<<" Puntos en la distro"<<endl;
  }

  cout<<"vamos bien."<<endl;

  PuntosY=PuntosX;
  int kuenta=0;

  for(int j=0; j<Puntos1D;j++){
    for(int i=j; i<Puntos1D;i++){   
      PuntosQ(kuenta,0)=PuntosX(j);
      PuntosQ(kuenta,1)=PuntosY(i);
      cout<<"Llevo "<<kuenta<<" Puntos en la distro completa"<<endl;
      kuenta++;
    }
  }
  
  PuntosQ.save("MuestreoPuntosOA2D.dat", raw_ascii);
 
  cout<<"Terminando de poblar boltzmanianamente el espacio "<<endl;
  double aux;

  for(int i=0; i<PuntosMalla-1;i++){
    for(int j=i+1; j<PuntosMalla;j++){  
      aux=UNormPropagatorFreeImagTime(PuntosQ(i,0), PuntosQ(j,0),dt)*
	UNormPropagatorFreeImagTime(PuntosQ(i,1), PuntosQ(j,1),dt);
      if(aux>0.00000000000001)propagator(i,j)=aux;
      propagator(j,i)=propagator(i,j);
    }
  }

  propagator=propagator/(dt*PuntosMalla);

  //  cout<<propagator;

  
  propagator.save("UNormPropagador02D.dat", raw_ascii);

  mat eigenestados;
  vec energias;
  
  eig_sym(energias, eigenestados, propagator);

  energias.save("LamdaExpEnergias.dat", raw_ascii);
  
  //energias=energias/(PuntosMalla*dt);

  energias=log(energias);

  energias=energias/(-dt);
  
  energias.save("Energias.dat", raw_ascii);
  
  eigenestados.save("EigenStates.dat", raw_ascii);

  return 0;

}

