//Version 1D para probar las ideas de eduardo.


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

  int Puntos1D=100;
  
  //double dt=1.0;
  if(argc==1){ 
    dt=0.1;
  }else{    
    dt=strtod(argv[1], NULL);
    Puntos1D=(int)strtod(argv[2],NULL);
    cout<<dt<<endl;
  }
  
  //  const double Z_V=pi*hbar*sqrt(2./omega)/dt; 

  

    cout<<"Probando exactamente lo de Eduardo, Empecemos con dos OA differentes."<<endl;

  double x,y;
  int contando=0; 
  
  vec PuntosX(Puntos1D);
  mat propagator=eye<mat>(Puntos1D,Puntos1D);

  const double xmin=-12, xmax=12;
  
  for(int i=0; i<Puntos1D;i++){
    //exactamente igual que en los adendums del paper
    //muestreo "perfecto"
    PuntosX(i)=xmin+((double)i+0.5)*(xmax-xmin)/(double)Puntos1D;
    
  }

  cout<<"Terminando de poblar homogeneamente el espacio "<<endl;
  double aux;
  rowvec quno,qdos;
  
  cout<<"Sacando el Propagador"<<endl;
  for(int i=0; i<Puntos1D;i++){
    quno=PuntosX.row(i);
    for(int j=0; j<Puntos1D;j++){
      qdos=PuntosX.row(j);
      
      //Eduardo Zambrano Ordering VTV
      aux=UNormPropagatorFreeImagTime(quno, qdos,dt)*
	PropOA1D(PuntosX(i), PuntosX(j),dt);
      if(aux>0.00001)propagator(i,j)=aux;     
      }
  }

  double area;
  area=(xmax-xmin);
  //propagator=propagator*(1.0/(PuntosMalla*2.*pi*hbar*dt));
  propagator=propagator*(1.0/sqrt(2.*pi*hbar*dt));
  propagator=propagator*(area/Puntos1D);
  

  propagator.save("PropEduardo1D.dat", raw_ascii);

   mat eigenestados;
  vec energias;
  
  cout<<"Haciendo la parte perra, sacando EigenEnergias"<<endl;
  eig_sym(energias, eigenestados, propagator);

  energias.save("LamdaExp1D.dat", raw_ascii);
  
  //energias=energias/(PuntosMalla*dt);

  energias=log(energias);

  energias=energias*hbar/(-dt);
  
  energias.save("Energias1D.dat", raw_ascii);
  
  eigenestados.save("EigenStates1D.dat", raw_ascii);

 
  return 0;

}

