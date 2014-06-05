/* Archivo que SOLO tiene constantes Globales para Mecanica Cuantica */
#ifndef __Q_CONST__
#define __Q_CONST__

const double pi=3.1415926;
//const double hbar=0.05;
const double hbar=1.0; //reproduciendo el paper
const double Energy=0.813840071;
/* dt tiene que respetar la desigualdad que esta en la segunda pagina del paper.
   1>>Z_V*(2 pi hbar)((N dt) 
   Resulta que Z_V es calculableen este caso y entonces dt approx 1 cumple la desigualdad , pero tiene que ser menor a 1 por el propagador que le falta O(dt^²)*/
//const double dt=0.2; 
const double dt=0.00549546;  //reproduciendo el paper
const double omega=0.1;
const double Z_V=pi*hbar*sqrt(2./omega)/dt; 
#endif
