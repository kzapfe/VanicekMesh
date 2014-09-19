/*
Vamos a intentar usar armadillo para obtener lo mismo que con 
Python, usando armadillo. Hold your breath, men.
 */

#include <armadillo> //intentemos primero con arma
#include <string>

using namespace std;
using namespace arma;


double PotNelson(double qx, double qy){
  double w1=0.1; 
  double result;
  result=(w1*qx*qx)/2.+(qy-qx*qx/2.)*(qy-qx*qx/2.);
  return result;

}


int main(int argc, char *argv[]){

  double dt=atof(argv[1]);
  string tailnombre;
  tailnombre=argv[1];

  const double hbar=0.05;
  const int N=256;
  const double xmin=-5.0, xmax=5.0;
  const double ymin=-1.0, ymax=9.0;
  const double pi=3.14159265359;
  

  string nombreenergias;
  string nombreestados;
  string colapuntos;
  stringstream ss;
  ss << N; 
  colapuntos = ss.str();

  nombreenergias="Energias-"+tailnombre+"-"+colapuntos+".dat";
  nombreestados="Estados-"+tailnombre+"-"+colapuntos+".dat";

  cout<<"Vamos a obtener los archivos "<<nombreestados<<
    " y "<<nombreenergias<<endl;


  //ugly but readable beats correct but illegible
  mat PuntosQ(N*N,2);
  
  int cuenta=0;
  
  for(int j=0; j<N; j++){
    for(int k=0; k<N; k++){
      PuntosQ(cuenta,0)=xmin+(xmax-xmin)*((double)j+1./2.)/(double)N;
      PuntosQ(cuenta,1)=ymin+(ymax-ymin)*((double)k+1./2.)/(double)N;
      cuenta++;
    }
  }

  //Hasta ahora, igualito que Python.
  PuntosQ.save("Coordenadas80.dat", raw_ascii);
  
  sp_mat Propagator(N*N, N*N);
  double valor=0.000;
  const double tolerancia=0.0001;
  rowvec Distancia(2);
  double auxiliar;

  for(int j=0; j<N*N; j++){
    for(int k=0; k<N*N; k++){
     
      Distancia=PuntosQ.row(j)-PuntosQ.row(k);
     
      auxiliar=dot(Distancia,Distancia);
     
      valor=exp(
		 -(auxiliar)/(2*dt*hbar)
		 -(PotNelson(PuntosQ(j,0),PuntosQ(j,1))+
		   PotNelson(PuntosQ(k,0),PuntosQ(k,1)))
		 *dt/(hbar*2.0)
		);
      
      if(abs(valor)>tolerancia){
	Propagator(j,k)=valor;
      }
      
    }
  }

  //Normalizamos el propagador
  double ZD=(xmax-xmin)*(ymax-ymin);
  Propagator=Propagator*ZD/(N*N*2.*hbar*dt*pi);

  int kuantos=1000;
  vec EigenValores(kuantos);
  mat EigenEstados;
  
  eigs_sym(EigenValores, EigenEstados, Propagator, kuantos);
  
  EigenValores=log(EigenValores)*hbar/(-dt);

  //voltear especularmente la matriz de estados  y las energias
  //para seguir con el estandar simple (energia baja, columna inmediata)
  EigenValores=flipud(EigenValores);
  EigenEstados=fliplr(EigenEstados);
  
  EigenEstados.insert_cols(0, PuntosQ);

  EigenValores.save(nombreenergias,raw_ascii);
  EigenEstados.save(nombreestados,raw_ascii);


  return 0;


}
