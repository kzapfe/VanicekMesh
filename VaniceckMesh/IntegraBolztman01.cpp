#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <gsl/gsl_rng.h> //Inicializa el gsl random
#include <gsl/gsl_randist.h> //las distribuciones de gsl


#include "QuantumConstants.hpp"
#include "BoltzmanforNelson.hpp"
#include "PropagatorFreeParticle.hpp"


int main(int argc, char *argv[]){


 
  const int PuntosMalla=100000;
 
  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  ofstream pancha,lola;
  pancha.open("PuntosyBoltzman.dat");

  double x,y;
  double NPot;
  double dado;
  double PDFBoltz;
  double integral=0;
  
  
  for(int i=0;i<PuntosMalla;i++){
    x=gsl_ran_flat( r, -7, 7);
    y=gsl_ran_flat( r, -2,14);  
    PDFBoltz=Boltzmanpdf(x,y);
    pancha<<x<<"\t"<<y<<"\t"<<PDFBoltz<<endl;
    integral+=PDFBoltz;    
    };

  integral=integral/PuntosMalla*14.*16.;

  cout<<"la integral de pdfBoltz al azar es "<<integral<<endl;
  
  integral=0.00;

  lola.open("PuntosIntegraBoltz01.dat");
  
  for(int i=0;i<500;i++){
    for(int j=0;j<500;j++){
      x=-7.+(double)i/500.*14.0;
      y=-2.+(double)j/500.*16.0;
    PDFBoltz=Boltzmanpdf(x,y);
    lola<<x<<"\t"<<y<<"\t"<<PDFBoltz<<endl;
    integral+=PDFBoltz;    
    };
  }
  
  integral=integral/500./500*14.*16.;
  
  cout<<"la integral de pdfBoltz al riemancuadros es "<<integral<<endl;


  return 0;

}
