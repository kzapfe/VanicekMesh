//Metodo de 2013 Vanicek et al para hacer super buenas 
//aproximaciones a estados cuanticos.
//First Test

//primero, todos los include no hechos por mi
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <gsl/gsl_rng.h> //Inicializa el gsl random
#include <gsl/gsl_randist.h> //las distribuciones de gsl
//#include <armadillo>


//#include <ctime> //para ver cuanto se tarda en cada punto
//#include <omp.h> //paralel

// despues todos los include de mi creacion

#include "QuantumConstants.hpp"
#include "BoltzmanforNelson.hpp"


using namespace std;

int main(int argc, char *argv[]){

  const int PuntosMalla=50000;
  const int NivelMaximo=100;
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);


  cout<<"Hola"<<endl;

  double x,y;
  double NPot;
  double dado;
  double PDFBoltz;
  int contando=0;

  while(contando<PuntosMalla){
    x=gsl_ran_flat( r, -7, 7);
    y=gsl_ran_flat( r, -2,14);
    dado=gsl_ran_flat( r, 0, 1.0);
    PDFBoltz=Boltzmanpdf(x,y);
    if(PDFBoltz>dado){
      NPot=NelsonPotential(x,y);      
      cout<<x<<"\t"<<y<<"\t"<<NPot<<"\t"<<PDFBoltz<<
	endl;
      contando++;
    }
  };

  

  return 0;

}
