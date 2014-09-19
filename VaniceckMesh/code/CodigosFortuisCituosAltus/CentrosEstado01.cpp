/* Programa que carga un estado de Nelson y saca los Centros en 
   el corte de las y=(qy,py) */

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
  const double xmin=-5.01,xmax=5.01;
  const double ymin=-1.01, ymax=9.01;
  const double pmin=-1.21, pmax=1.21;

  mat DominioQ(PuntosDominioEje,2);
  double dx=(xmax-xmin)/PuntosDominioEje;
  double dy=(ymax-ymin)/PuntosDominioEje;
  double dp=(pmax-pmin)/PuntosDominioEje;
  
  //El dominio en x no te importa una vez integrado
  for(int j=0; j<PuntosDominioEje; j++){
    DominioQ(j,0)=ymin+j*dy;
    DominioQ(j,1)=pmin+j*dp;
  }

  // DominioQ.save("DominioQ.dat", raw_ascii);

  /*Now comes the fun part */

  /*select section constant x */
  const double qx=0, px=0;
  vec Integradoxhiqx=zeros<vec>(PuntosDominioEje);
  int auxiliar;

  for(int j=0; j<PuntosDominioEje; j++){
    Integradoxhiqx(j)=dot(valores.row(j),fliplr(valores.row(j)));
  }

  mat test(PuntosDominioEje,2);
  
  test.col(0)=DominioQ.col(0);
  test.col(1)=Integradoxhiqx;
  
  test.save("TestingAutoCorrelations.dat", raw_ascii);

  double centroq;
  double centrop;
  double cuerdaaux;
  double integrando=0;
  mat Wigner(PuntosDominioEje,PuntosDominioEje);
  Wigner.zeros();
  int limiteintegracuerdas;
  double exponente=0;

  /* A los indices de abajo les quitaste el primer y ultimo
     puntos, ya que no podemos resolver la integral hasta la orilla, 
     sino solo hasta un punto antes */
  for(int j=1; j<PuntosDominioEje-1; j++){
    //indice de las qy
    centroq=ymin+j*dy;
    limiteintegracuerdas=min(j,PuntosDominioEje-1-j);

    for(int k=1; k<PuntosDominioEje-1; k++){
      //indice de las py
      centrop=pmin+k*dp;
      
      for(int l=-limiteintegracuerdas; l<=limiteintegracuerdas; l++){
	//integramos sobre las cuerdas
	/*Aqui es donde realmente calculamos la transformada para Centros */

	cuerdaaux=l*dy;
	//exponente=l*k      ponente=(-2.*centrop*cuerdaaux/hbar);
	//columna j-> Posicion q, 
	//renglon j -> Momento p
	Wigner(k,j)+=Integradoxhiqx(j-l)*Integradoxhiqx(j+l)*cos(exponente);
	//cout<<j<<"\t"<<k<<"\t"<<l<<"\t"<<cos(exponente)<<"\t"<<endl;
      }
      //cout<<j<<"\t"<<k<<"\t"<<"\t"<<Wigner(k,j)<<"\t"<<endl;     
    }
  }
     
  //las dy tienen la misma escala que las d(cuerdaauxiliar)
  Wigner=Wigner*dx*dy/(sqrt(2.*pi*hbar));

  Wigner.save("TestWigner.dat",raw_ascii);

  return 0;
  
};

