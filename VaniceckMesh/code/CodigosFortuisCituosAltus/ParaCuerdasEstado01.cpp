/* Programa que carga un estado de Nelson y saca las Cuerdas  en 
   el corte conjugado a  las y=(qy,py) (xhi_2)*/
/* Version en Paralelo*/

const double hbar=0.05;
const double pi=3.1415926;

#include<cmath>
#include<armadillo>

using namespace std;
using namespace arma;


int main(int argc, char *argv[]){

  if ( argc != 2 ) 
    {// argc should be 2 for correct execution
    // We print argv[0] assuming it is the program name
      cout<<"Uso: "<< argv[0] <<" <Nombre_de_archivo_estado_psi(q)>\n";
      exit(1);
    }
  
  string NombreValores;
  NombreValores=argv[1];
  mat valores;
  valores.load(NombreValores);

  const int PuntosDominioEje=valores.n_rows;
  cout<<PuntosDominioEje<<endl;

  const double xmin=-5.01,xmax=5.01;
  const double ymin=-1.01, ymax=9.01;
  const double pmin=-1.21, pmax=1.21;
  
  //la mitad del dominio en muq es  la mitad del TAMANHO del
  //dominio de las yq, istus est:
  const double muqmin=-5.01;
 
  mat DominioQ(PuntosDominioEje,2);
  double dx=(xmax-xmin)/PuntosDominioEje;
  double dy=(ymax-ymin)/PuntosDominioEje;
  double dp=(pmax-pmin)/PuntosDominioEje;
  
 
  
  /*
  //El dominio en x no te importa una vez integrado
  for(int j=0; j<PuntosDominioEje; j++){
    DominioQ(j,0)=ymin+j*dy;
    DominioQ(j,1)=pmin+j*dp;
  }
  */
  
  const double dmuq=dy;
  const double dmup=dp;


  //Date cuenta que tus cuerdas no pueden tener
  //mas resolucion que tus centros. De hecho, tienen
  //la mitad de resolucion, por la definicion 
  //de tus centros. Tal vez puedas interpolar mas finamente,
  //pero recuerda: al final siguen siendo 
  //putos 100 por 100 puntos originalmente!!!
 
  

  ofstream cuerdasout;
  cuerdasout.open("CuerdasEnMu.dat");

  /* A los indices de abajo les quitaste el primer y ultimo
     puntos, ya que no podemos resolver la integral hasta la orilla, 
     sino solo hasta un punto antes */
  
  /* Escogamos la cuerda mu, Primero*/

  int PuntosDominioMu=PuntosDominioEje/8-1;
  //De momento, vamos a centrarnos cerca del origen, asi que no vamos
  // a calcular todo. TODO es PuntosDominioEje/2, por simetria.
  //si queremos un tercio... sera entre 6. No podemos tener "mas informacion".

  mat CuerdaEnDosXhiReal=zeros<mat>(2*PuntosDominioMu, 2*PuntosDominioMu);
  mat CuerdaEnDosXhiImag=zeros<mat>(2*PuntosDominioMu, 2*PuntosDominioMu);

  cout<<PuntosDominioMu<<endl;

 
#pragma omp parallel num_threads(8)
  {
    

#pragma omp  for schedule(dynamic)
    //Se supone que todo lo que esta declarado dentro del for es privado.
    
    //es mas simple si usamos una notacion simetrica, aunque no sea lo 
    //mas c++esco. Readibility first
    for(int j=-PuntosDominioMu; j<PuntosDominioMu; j++){
      //j va a contar los indices en mup
      double muq;
      double mup;
      double yq;
      double integrando=0;
      int limiteintegrayq;
      double exponente=0;   

      mup=j*dmup;
    
      for(int k=-PuntosDominioMu; k<PuntosDominioMu; k++){
      //indice de las muq 
      //muq determina hasta donde vamos a poder integrar en q_y
	muq=k*dmuq;
     
	limiteintegrayq=PuntosDominioEje-abs(k)-1;
	
	double auxiliarresult=0;
	
	for(int l=abs(k); l<=limiteintegrayq; l++){
	  //en este ciclo integramos sobre y.
	  yq=ymin+l*dy;
	  //esto efectivamente integra sobre las x:
	  //cout<<"probable lugar de cagarla"<<endl;
	  //cout<<l<<"\t"<<l+k<<"\t"<<l-k<<"\t"<<k<<endl;
	  
	  auxiliarresult=dot(valores.row(l-k),valores.row(l+k));
	  exponente=(2.*yq*mup/hbar);
	//selecciona el renglon y-muq y el y+muq:
	  CuerdaEnDosXhiReal.at(k+PuntosDominioMu,j+PuntosDominioMu)
	    +=auxiliarresult*cos(exponente);
	  CuerdaEnDosXhiImag.at(k+PuntosDominioMu,j+PuntosDominioMu)
	    +=auxiliarresult*sin(exponente);
	  
	}//se termina la integracion numerica
	
     
     
      }
    
    
    }
    
  } //close omp pragma parallel 

  double muq;
  double mup;



   CuerdaEnDosXhiImag*=(dy*dx)/(2.*pi*hbar);
   CuerdaEnDosXhiReal*=(dy*dx)/(2.*pi*hbar);

   for(int j=-PuntosDominioMu; j<PuntosDominioMu; j++){
     //j va a contar los indices en mup
    mup=j*dmup;
    
    for(int k=-PuntosDominioMu; k<PuntosDominioMu; k++){
	//indice de las muq 
      //muq determina hasta donde vamos a poder integrar en q_y
      muq=k*dmuq;
  
      cuerdasout<<2*muq<<"\t"<<2*mup<<"\t"
		<<CuerdaEnDosXhiReal(k+PuntosDominioMu,j+PuntosDominioMu)
		<<"\t"<<CuerdaEnDosXhiImag(k+PuntosDominioMu,j+PuntosDominioMu)
		<<endl;     

    }
    cuerdasout<<endl;
   }


  return 0;
  
};

