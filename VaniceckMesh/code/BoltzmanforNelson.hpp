//la distribucion de Boltzman para el potencial de Nelson
#ifndef __BoltzmanForNelson__
#define __BoltzmanForNelson__


using namespace std;

double PotArm(double x, double y){
  double result;
  //Caso uno, isotropico, degenerado. 
  const double w1=1.0, w2=1.0;
  result=w1*x/2.+w2*y/2;

  return result;

}

double NelsonPotential(double x, double y){
  //parametros de Nelson 
  const double omega1=0.1;
  const double omega2=1;  
  
  double result;

  result=omega1*x*x/2.0+omega2*(y-x*x/2.)*(y-x*x/2.);
  
  return result;



};

double Boltzmanpdf(double x, double y){

  double result;
  result=exp(-dt*NelsonPotential(x,y)/hbar);  
  return result;

}



#endif

