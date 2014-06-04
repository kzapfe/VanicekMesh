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
    x=gsl_ran_flat( r, -14, 14);
    y=gsl_ran_flat( r, -6,18);  
    PDFBoltz=Boltzmanpdf(x,y);
    pancha<<x<<"\t"<<y<<"\t"<<PDFBoltz<<endl;
    integral+=PDFBoltz;    
    };

  integral=integral/PuntosMalla*28.*24.;

  cout<<"la integral de pdfBoltz al azar es "<<integral<<endl;
  
  integral=0.00;

  lola.open("PuntosIntegraBoltz01.dat");
  
  for(int i=0;i<700;i++){
    for(int j=0;j<700;j++){
      x=-14.+(double)i/500.*28.0;
      y=-6.+(double)j/500.*24.0;
    PDFBoltz=Boltzmanpdf(x,y);
    lola<<x<<"\t"<<y<<"\t"<<PDFBoltz<<endl;
    integral+=PDFBoltz;    
    };
  }
  
  integral=integral/500./500*28.*24.;
  
  cout<<"la integral de pdfBoltz al riemancuadros es "<<integral<<endl;
  cout<<"El valor analitico deberia dde ser  "<<Z_V<<endl;


  return 0;

}
