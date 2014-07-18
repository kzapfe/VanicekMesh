//Metodo de 2013 Vanicek et al para hacer super buenas 
//aproximaciones a estados cuanticos.
//First Test

//primero, todos los include no hechos por mi
#include <iostream>
#include <fstream>
#include <vector>
//#include <armadillo>

//#include <ctime> //para ver cuanto se tarda en cada punto
#include <omp.h> //para ver cuanto se tarda en cada punto

// despues todos los include de mi creacion

#include "QuantumConstants.hpp"
#include "BoltzmanforNelson.hpp"


using namespace std;

int main(int argc, char *argv[]){

  const int PuntosMalla=500;
  const int NivelMaximo=100;

  cout<<"Hola"<<endl;
  cout<<"valores de Nelson a lo guey"<<endl;
  
  cout<<NelsonPotential(0,1)<<endl;
  cout<<NelsonPotential(2,1)<<endl;

  cout<<NelsonPotential(2,2)<<endl;
  cout<<NelsonPotential(9,3)<<endl;
  
  


  return 0;

}
