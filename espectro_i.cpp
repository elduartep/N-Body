//Calcula el espectro de potencias dadas las posiciones inciales de las particulas

//paquetes gnu
#include<iostream>
#include<fstream>
#include<cmath>
#include<ctime>
#include <alloca.h>
//FFTW
#include <fftw3.h>
#include <complex.h>
//numerical recipes
#include "nr3.h"
//mis bloques
#include "funciones.h"
#include "otro.h"
#include "crece.h"


using namespace std;


int main(){

////////////////////////////		FFTW en paralelo
int fftw_init_threads (void); 
cout<<"debe ser direfrente de cero "<<fftw_init_threads ()<<endl;
void fftw_plan_with_nthreads (int nthreads);
fftw_plan_with_nthreads (4);
// -lfftw3_threads -lfftw3 -lm -lpthread
// -lfftw3_omp -lfftw3 -lm -fopenmp

	clock_t t1,t2;
	t1=clock();

	int aux2,mue,nd;
	int i,j,k,j2,k2,l;
	int p, ii, jj, kk, iii, jjj, kkk, ijk; //los dobles y triples son i+-1, j+-1, k+-1.    ijk es un contador que barre el numero de celdas

	double a,h,cK;					
	double gxp, gyp, gzp;
	double dx, dy, dz, tx, ty, tz; //variables auxliares para calcular la contribucion de cada particula a las celdas adyacentes
	double aux, auxk;
	double elk, kr, r1, r2, kmin, kmax;
	double ff;//es la funcion f
	double a1,a2,a3,a4;
	double sx,sy,sz;//auxiliares




	mm  = (long int*) malloc(sizeof(long int) * bin);
	MpI = (double*) fftw_malloc(sizeof(double) * bin );
	MpR = (double*) fftw_malloc(sizeof(double) * bin );
	Mk  = (double*) fftw_malloc(sizeof(double) * bin );
	EpI = (double*) fftw_malloc(sizeof(double) * bin );
	EpR = (double*) fftw_malloc(sizeof(double) * bin );
	Ek  = (double*) fftw_malloc(sizeof(double) * bin );

	rho = (double*) fftw_malloc(sizeof(double) * n3 );
	delta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n2* (n/2+1));
	dft = fftw_plan_dft_r2c_3d(n, n, n, rho, delta, FFTW_MEASURE);

	x = (double*) fftw_malloc(sizeof(double) * np3 );
	y = (double*) fftw_malloc(sizeof(double) * np3 );
	z = (double*) fftw_malloc(sizeof(double) * np3 );

	char nomArch[22];
	FILE * NOM[22];
	FILE * O[n_mue*n_k];
	FILE * E;
	FILE * L;


double Di,Dpi;//la funcion de crecimiento y su derivada en (ai) y (ai+da/2) respectivamente
cout<<"funcion de crecimiento"<<endl;
crece(Di,Dpi);


	for(mue=mue_i;mue<mue_f;mue++){//	ciclo sobre el numero de muestras	
	cout<<"calculando el espectro conjunto para la muestra "<<muestra[mue]<<endl;
	aux2=0;//for(aux2=z_espectro;aux2<6;aux2++){	//	ciclo sobre z

	 	E=fopen("temp.dat","w+");


		for(l=k_i;l<k_f;l++){				//	ciclo sobre el numero de onda
			if(aux2==5){//es para calcular el espectro de los halos, -> solo para z=0
				sprintf(nomArch,"espectro%c_%d.dat",muestra[mue],l);
				O[mue*n_k+l]=fopen(nomArch,"w+");}



/////////////////////////////         leyendo el archivo con las posiciones y velocidades iniciales
	sprintf(nomArch,"in%c%d_%d.dat",muestra[mue],l,np);
	NOM[0]=fopen(nomArch,"r");
	ijk=0;
	aux=pow(ai-aa*0.5,2.)*Dpi*alpha[l];
	double con=1.*np/n+n*gama;
	double cons=Di*alpha[l];
	for(k=0;k<np;k++){
	  for(j=0;j<np;j++){
	    for(i=0;i<np;i++){
	    	fscanf(NOM[0],"%le %le %le\n",&sx,&sy,&sz);//escanea las posiciones iniciales
				//	asigna las posiciones iniciales usando la aproximacion de Zeldovich
		    x[ijk]=fmod(con+1.*i*n/np-cons*sx,n);
		    y[ijk]=fmod(con+1.*j*n/np-cons*sy,n);
		    z[ijk]=fmod(con+1.*k*n/np-cons*sz,n);
//if((abs(vx[ijk])>control) || (abs(vy[ijk])>control) || (abs(vz[ijk])>control)){
//cout<<endl<<a<<"	"<<ijk<<"	alerta,	delta_xyz=	"<<vx[ijk]<<"	"<<vy[ijk]<<"	"<<vz[ijk]<<"	"<<endl;
//cout<<"delta_a= "<<(abs(vx[ijk])+abs(vy[ijk])+abs(vz[ijk]))/control<<endl<<endl;}
				ijk++;}}}
  fclose(NOM[0]);







////////////////////////////////////////////		calculo de la densidad

			for(i=0;i<n3;i++)//loop sobre las celdas para borrar la densidad antigua
			  rho[i]=-1.;

			for(p=0;p<np3;p++){//loop sobre las particulas para calcular la densidad de las celdas
			  i=floor(x[p]);		j=floor(y[p]);		k=floor(z[p]);
			  dx=x[p]-1.*i;		dy=y[p]-1.*j;		dz=z[p]-1.*k;
			  tx=1.-dx;		ty=1.-dy;			tz=1.-dz;
			//condiciones de frontera periodicas
			  ii=(i+1)%n;		jj=(j+1)%n;		kk=(k+1)%n;	
			//incrementando la densidade de las celdas adyacentes	(solo las componentes reales)	
			  rho[i+j*n+k*n2]+=m*tx*ty*tz;		rho[ii+j*n+k*n2]+=m*dx*ty*tz;
			  rho[i+jj*n+k*n2]+=m*tx*dy*tz;		rho[ii+jj*n+k*n2]+=m*dx*dy*tz;
			  rho[i+j*n+kk*n2]+=m*tx*ty*dz;		rho[ii+j*n+kk*n2]+=m*dx*ty*dz;
			  rho[i+jj*n+kk*n2]+=m*tx*dy*dz;	rho[ii+jj*n+kk*n2]+=m*dx*dy*dz;}



		//este es el calculo de la transformada de Fourier de la sobredensidad
			fftw_execute(dft);


	
		//imprimiendo el espectro calculado
			for(k=0;k<n;k++){//loop sobre los numeros de onda
			  for(j=0;j<n;j++){
			    for(i=0;i<n/2+1;i++){
		  	    ijk=i+(n/2+1)*(j+k*n);
						if((j==0) || (j==n/2))	{j2=j;}	else	{j2=n-j;}		if(j<j2)	{j2=j;}		else	;
						if((k==0) || (k==n/2))	{k2=k;}	else	{k2=n-k;}		if(k<k2)	{k2=k;}		else	;
						elk=dpin*sqrt(k2*k2+j2*j2+i*i);
						//if((elk>0.1) && (elk<1.)){
						if(elk<corte){
							fprintf(E,"%le %le %le\n",elk*e[l],delta[ijk][0]*delta[ijk][0]/pow(e[l]*n,3.),delta[ijk][1]*delta[ijk][1]/pow(e[l]*n,3.));
							if(aux2==5)
								fprintf(O[mue*n_k+l],"%le %le %le\n",elk*e[l],delta[ijk][0]/pow(e[l]*n,1.5),delta[ijk][1]/pow(e[l]*n,1.5));}}}}

		if(aux2==5)
			fclose(O[mue*n_k+l]);

		}	//	cierra el ciclo sobre el numero de onda

	  fclose(E);








//////////////////////////////////		haciendo una media sobre los numeros de onda

//calcula la media
		L=fopen("temp.dat","r");
		nd=0;
		kmin=100.;
		kmax=-1.;
		while(fscanf(L,"%le	%le	%le\n",&kr, &r1, &r2)!=EOF){
			nd++;										//	cuenta el numero de lineas de entrada
			if((kmin>kr) && (kr!=0.))	kmin=kr;		//	busca el limite inferior de k
			if(kmax<kr)	kmax=kr;}		//	busca el limite superior de k

		FILE * Km;
		Km=fopen("Km.dat","w+");
		fprintf(Km,"%le	%le\n",kmin,kmax);
		fclose(Km);

		L=fopen("temp.dat","r");
		int entero=pow(10,8);
		cK=(log(entero*kmax)/log(entero*kmin)-1.)/bin;

//borra el contador, la media y el error
		for(i=0;i<bin;i++){ MpI[i]=0.;	MpR[i]=0.;	Mk[i]=0.;	EpI[i]=0.;	EpR[i]=0.;	Ek[i]=0.;	mm[i]=0;}

//incrementa usando los datos
		for(i=0;i<nd;i++){
			fscanf(L,"%le	%le	%le\n",&kr, &r1, &r2);
			for(j=0;j<bin;j++){
				if( (kr>=pow(entero*kmin,1.+cK*j)/entero) && (kr<pow(entero*kmin,1+cK*(j+1))/entero) ){
//				if((kr>(kmin*pow(kmax/kmin,j*1./bin))) && (kr<(kmin*pow(kmax/kmin,(j*1.+1.)/bin)))){
					mm[j]++;		MpR[j]+=r1;		MpI[j]+=r2;		Mk[j]+=kr;}}}

//normaliza dividiendo por el numero de ocurrencias
		for(i=0;i<bin;i++){MpR[i]/=mm[i];	MpI[i]/=mm[i];	Mk[i]/=mm[i];}

//calcula el desvio al cuadrado
		L=fopen("temp.dat","r");
		for(i=0;i<nd;i++){
			fscanf(L,"%le	%le	%le\n",&kr, &r1, &r2);
			for(j=0;j<bin;j++)
				if( (kr>=pow(entero*kmin,1.+cK*j)/entero) && (kr<pow(entero*kmin,1+cK*(j+1))/entero) ){
//				if((kr>(kmin*pow(kmax/kmin,j*1./bin))) && (kr<(kmin*pow(kmax/kmin,(j*1.+1.)/bin)))){
		  	  EpR[j]+=pow(r1-MpR[j],2);
		  	  EpI[j]+=pow(r2-MpI[j],2);
		  	  Ek[j]+=pow(kr-Mk[j],2);}}
		fclose(L);


//calcula sigma²
		for(i=0;i<bin;i++){
		  if(mm[j]==0 || mm[j]==1)
				cout<<"error, hay "<<mm[j]<<" halos en el bin "<<j<<endl;
	    else{
				EpR[i]/=(mm[i]*(mm[i]-1));
			  EpI[i]/=(mm[i]*(mm[i]-1));
			  Ek[i]/=(mm[i]*(mm[i]-1));}}


//sigma, y normalizando
		sprintf(nomArch,"espectro%c_%d_i.dat",muestra[mue],aux2);
		NOM[aux2+11]=fopen(nomArch,"w+");
		for(j=0;j<bin;j++){
	    aux=(MpI[j]+MpR[j]);
	    MpI[j]=aux;
	    aux=sqrt(EpI[j]+EpR[j]);
	    auxk=sqrt(Ek[j]);
	    fprintf(NOM[aux2+11],"%le %le %le %le\n",Mk[j],MpI[j],auxk,aux);}
		fclose(NOM[aux2+11]);



//	}		//	cierra el ciclo sobre z
	}		//	cierra el ciclo de as muestras	

	fftw_free(z);		fftw_free(y);		fftw_free(x);						//posicion de las particulas
	fftw_destroy_plan(dft);
	fftw_free(delta);	fftw_free(rho);

	fftw_free(MpI);		fftw_free(MpR);		fftw_free(Mk);		fftw_free(EpI);
	fftw_free(EpR);		fftw_free(Ek);		free(mm);











	cout<<"calculando la media de las muestras"<<endl;
//////////////////////////////////		haciendo una media de las muestras
	Ek = (double*) fftw_malloc(sizeof(double) * bin );
	EpI = (double*) fftw_malloc(sizeof(double) * bin * 2 );
	Mk = (double*) fftw_malloc(sizeof(double) * bin );
	MpI = (double*) fftw_malloc(sizeof(double) * bin );


aux2=0;//	for(aux2=z_espectro;aux2<6;aux2++){							//	ciclo sobre z
		for(i=0;i<bin;i++){ MpI[i]=0.;	Mk[i]=0.;	EpI[i]=0.;	EpI[bin+i]=0.;	Ek[i]=0.;}	//	borra las medias
		//hace la media
		for(mue=mue_i;mue<mue_f;mue++){		//	ciclo sobre las muestras
			sprintf(nomArch,"espectro%c_%d_i.dat",muestra[mue],aux2);
			NOM[aux2]=fopen(nomArch,"r");
			for(j=0;j<bin;j++){
		    fscanf(NOM[aux2],"%le %le %le %le\n",&a1,&a2,&a3,&a4);
				Mk[j]=a1;			MpI[j]+=a2;			Ek[j]=a3;			EpI[j]+=a4*a4;}
			fclose(NOM[aux2]);}
		for(j=0;j<bin;j++)
			MpI[j]/=n_mue;

		//	calcula el desvio al cuadrado
		for(mue=mue_i;mue<mue_f;mue++){		//	ciclo sobre las muestras
			sprintf(nomArch,"espectro%c_%d_i.dat",muestra[mue],aux2);
			NOM[aux2]=fopen(nomArch,"r");
			for(j=0;j<bin;j++){
		    fscanf(NOM[aux2],"%le %le %le %le\n",&a1,&a2,&a3,&a4);
	  	  EpI[bin+j]+=pow(a2-MpI[j],2);}
			fclose(NOM[aux2]);}

		//calcula sigma²
		for(i=0;i<bin;i++){
		  EpI[bin+i]/=(n_mue*n_mue);
		  EpI[i]/=n_mue;
			EpI[i]+=EpI[bin+i];}


//sigma
		sprintf(nomArch,"espectro_%d_i.dat",aux2);
		NOM[aux2+11]=fopen(nomArch,"w+");
		for(j=0;j<bin;j++){
	    aux=sqrt(EpI[j]);
	    fprintf(NOM[aux2+11],"%le %le %le %le\n",Mk[j],MpI[j],Ek[j],aux);}
		fclose(NOM[aux2+11]);
//	}		//	cierra el ciclo sobre z

cout<<"antes de cerrar"<<endl;
			

	fftw_free(MpI);		fftw_free(Mk);		fftw_free(EpI);		fftw_free(Ek);



	cout<<"acaba espectro"<<endl<<endl;

t2=clock();
double diff =((double)t2-(double)t1)/CLOCKS_PER_SEC;
cout<<diff<<" segundos,	"<<diff/60<<" minutos"<<endl;



}

