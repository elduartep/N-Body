//	Particle Mesh for N-Body evolution

//gnu
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <alloca.h>
//OPENMP
#include <omp.h>	
//FFTW
#include "fftw3.h"
//numerical recipes
#include "nr3.h"
#include "ran.h"
#include "gamma.h"
#include "deviates.h"
#include "interp_1d.h"
//mis bloques
#include "parameters.h"
#include "funciones.h"
#include "otro.h"
#include "inicio.h"
#include "crece.h"
#include "anda.h"
//#include "pelicula.h"		// modified anda.h to print the density field of the z=0 sheet 

using namespace std;

int main(){

cout<<endl<<"numer of particles = "<<np<<"^3"<<endl;
cout<<"numero de cells for Particle Mesh = "<<n<<"^3"<<endl;
cout<<"Lbox = "<<Lbox<<" Mpc/h"<<endl;


////////////////////////////		FFTW en paralelo
int fftw_init_threads (void); 
cout<<fftw_init_threads ()<<" must be different from zero "<<endl;
void fftw_plan_with_nthreads (int nthreads);
fftw_plan_with_nthreads (nt);
// -lfftw3_threads -lfftw3 -lm -lpthread
// -lfftw3_omp -lfftw3 -lm -fopenmp
omp_set_num_threads(nt);
cout<<"number of threads to be used = "<<nt<<endl;


int i;

//	PRIMERO LLAMAMOS A INICIO (internamente el barre las escalas y el numero de muestras indicadas)

inicio();

//	CALCULA LA FUNCION DE CRECIMIENTO
double D_i,Dp_i;//la funcion de crecimiento y su derivada en (ai)
double D_antes,Dp_antes;//la funcion de crecimiento y su derivada en(ai-da/2) respectivamente
crece(D_i,Dp_i,D_antes,Dp_antes);

//	AHORA QUEREMOS EVOLUCIONAR
anda(D_i,Dp_i,D_antes,Dp_antes);	//archivo	Di	Dpi



return 0;

}




//	compute the power spectrum by using Bias_pragma_lb.c, input 512x256_03.txt


