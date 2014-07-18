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


#include <ctime> //para ver cuanto se tarda en cada punto
//#include <omp.h> //paralel

// despues todos los include de mi creacion

#include "QuantumConstants.hpp"
#include "BoltzmanforNelson.hpp"
#include "PropagatorFreeParticle.hpp"



using namespace std;
using namespace arma;

int main(int argc, char *argv[]){

  int Puntos1D=20;
  
  //double dt=1.0;
  if(argc==1){ 
    dt=0.1;
  }else{    
    dt=strtod(argv[1], NULL);
    Puntos1D=(int)strtod(argv[2],NULL);
    cout<<dt<<endl;
  }
  
  //  const double Z_V=pi*hbar*sqrt(2./omega)/dt; 

  const int PuntosMalla=Puntos1D*(Puntos1D);
  //En este test, vamos a usar hbar=1

  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);


  cout<<"Probando exactamente lo de Eduardo, Empecemos con dos OA differentes."<<endl;

  double x,y;
  int contando=0; 
  mat propagator=zeros<mat>(PuntosMalla,PuntosMalla);
  vec PuntosX(Puntos1D), PuntosY(Puntos1D);
  mat PuntosQ(PuntosMalla,2);

  const double xmin=-4, xmax=4;


  for(int i=0; i<Puntos1D;i++){
    //exactamente igual que en los adendums del paper
    //muestreo "perfecto"
    PuntosX(i)=xmin+((double)i+0.5)*(xmax-xmin)/(double)Puntos1D;
    
  }

  cout<<"vamos bien."<<endl;

  PuntosY=PuntosX;
  double ZD=(xmax-xmin)*(xmax-xmin);
  

  int kuenta=0;

  for(int j=0; j<Puntos1D;j++){
    for(int i=0; i<Puntos1D;i++){   
      PuntosQ(kuenta,0)=PuntosX(j);
      PuntosQ(kuenta,1)=PuntosY(i);
      //  cout<<"Llevo "<<kuenta<<" Puntos en la distro completa"<<endl;
      kuenta++;
    }
  }
  
  PuntosQ.save("MuestreoPuntosOA2D.dat", raw_ascii);
  cout<<"Terminando de poblar homogeneamente el espacio "<<endl;
  double aux;
  rowvec quno,qdos;

  clock_t start, finish;
  
  start=clock();
  cout<<"Sacando el Propagador"<<endl;
  for(int i=0; i<PuntosMalla-1;i++){
    quno=PuntosQ.row(i);
    for(int j=0; j<PuntosMalla;j++){
      qdos=PuntosQ.row(j);
      
      aux=UNormPropagatorFreeImagTime(quno, qdos,dt)
	*UNormPropOAITime(PuntosQ(i,0), PuntosQ(j,0), 
			  PuntosQ(i,1), PuntosQ(j,1), dt) ;
      if(aux>0.00001)propagator(i,j)=aux;
    }
  }
  finish=clock();
  cout<<"Nos tardamos  "<< finish-start 
      << " ciclos de cpu en obtener el propagador"<<endl;;
 

   //propagator=propagator*(1.0/(PuntosMalla*2.*pi*hbar*dt));
  propagator=propagator*(1.0/(2.*pi*hbar*dt));
  propagator=propagator*(ZD/PuntosMalla);
 

  //  cout<<propagator;
  
  //  propagator.save("PropEduardoOA2D.dat", raw_ascii);

  mat eigenestados;
  vec energias;
  
  cout<<"Haciendo la parte perra, sacando EigenEnergias"<<endl;
  start=clock();
  eig_sym(energias, eigenestados, propagator);
   finish=clock();
   cout<<"Nos tardamos  "<< finish-start 
       << " ciclos de cpu en obtener la EigenDescomposicion"<<endl;;
   

  energias.save("LamdaExpEnergias.dat", raw_ascii);
  
  //energias=energias/(PuntosMalla*dt);

  energias=log(energias);

  energias=energias*hbar/(-dt);
  
  energias.save("Energias.dat", raw_ascii);
  
  eigenestados.save("EigenStates.dat", raw_ascii);

  cout<<"El estado base tiene energia " << energias(PuntosMalla-1)<<endl;

  return 0;

}

