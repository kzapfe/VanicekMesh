/* programa que suma varios estados coherentes sobre un oscilador
   armonico, representacion de Wigner */


#include <armadillo>
#include <iostream>
#include <cstdlib>

using namespace arma;
using namespace std;

using std::cout;


mat Rt(double w, double t){
  //Matriz de Rotacion aplastada
  mat Rt(2,2);
    Rt<<cos(w*t)<<-w*sin(w*t)<<endr
    <<sin(w*t)/w<<cos(w*t)<<endr;
    return Rt;
};


mat Mt(double w, double t){
  //Matriz de monodromía
  mat Mt(2,2);
    Mt<<1.0<<0.0<<endr
    <<0.0<<1.0<<endr;
    Mt=Mt+strans(Rt(w,t))*Rt(w,t);
    return Mt;
}; 


int main(){
  

  // Constantes
  int const N=3;
  double const chbarra=0.5;
  mat Jota;
  
  Jota<<0.0<<-1.0<<endr
		<<1.0<<00<<endr;

  

  //Familia de vectores sobre la cual hacer la suma
  vec eta[N]=ones<vec>(2);
  vec sigma[N][N];
  vec delta[N][N];
  
  eta[0]<<0.0<<endr
	<<0.0<<endr;
  eta[1]<<0.0<<endr
	<<1.0<<endr;
  eta[2]<<1.0<<endr
	<<0.0<<endr;
  
  
  complex<double> Ncuad=0.00;
  complex<double> esponente=0.00;
  

  //Constante de 
  for(int i=0; i<N;i++){
    for(int j=0; j<N;j++){
      esponente.real()=-norm(eta[i]-eta[j],2)/chbarra;
      esponente.imag()=dot(Jota*eta[i],eta[j])/chbarra;
      Ncuad=+exp(esponente);
    };};
  Ncuad=pow(Ncuad,-2);
      


  //inicializar los vectores suma y diferencia
  for(int i=0; i<N;i++){
    for(int j=0; j<N; j++){
      sigma[i][j]=eta[i]+eta[j];
      delta[i][j]=eta[i]-eta[j];
    };
  };
  
  /*La funcion F(t,w) solo nos interesa para
    argumentos reales positivos*/


 for(double w=0.011;w<3.0;w+=0.01){
  
 for(double t=0.00;t<8.10;t+=0.01){
      complex<double> Func=0.0;
      //Los beta solo aparecen como el argumento de la exp
      //Asi que no requieren indices
      cx_vec beta(2);


      //iniciar el vector beta
      for(int i=0; i<N;i++){
	for(int j=0; j<N; j++){
	  for(int k=0; k<N; k++){
	    for(int n=0; n<N; n++){
	      //beta solo aparece como el preliminario de Func
	      // asi que puedes tirar a la basura los indices
	      beta.set_real(delta[i][j]-strans(Rt(w,t))*delta[k][n]);	
	      beta.set_imag(Jota*sigma[i][j]-strans(Rt(w,t))*Jota*sigma[k][n]);	      
	      Func=+exp(as_scalar(strans(beta)*inv(Mt(w,t))*beta)/chbarra);	      
	    };};};};
    
      Func=2.0*Func*(double)pow(N,4)/sqrt(det(Mt(w,t)));

	cout<<t<<"\t"<<w<<"\t"<<Func.real()<<"\t"<<Func.imag()<<endl;

    };
    cout<<endl;
    };

      
  

  return 0;
    
}



