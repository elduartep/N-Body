
const double pi=4.*atan(1.);				//	numero pi
const double pi_2=1/(4.*atan(1.)*4.*atan(1.));		//	pi^-2


double dpin=2*pi/n;
double dpinp=2*pi/Lbox;//*Lbox/n;
int np3=np*np*np;
int np2=np*np;
int n3=n*n*n;
int n2=n*n;


////////////////////////////		main variables
	double *rho;					//	density contrast
	double *x, *y, *z;				//	particles' position
	double *vx, *vy, *vz;				//	particles' velocity
	double *gx, *gy, *gz;				//	acceleration in each cells center
	double *MpI, *MpR, *Mk, *EpI, *EpR, *Ek;	//	variables for computing P(k)
	double *Mkb, *Mpb, *Ekb, *Epb;
	long int *mm;
	double *Mmin_k,*Mmax_k;				//	max and min mass for each box size
	double *vc;					//	correlation fuction
	double *rc;					//	radius for the correlation function
	long int *ic;					//	numbercount




//	FUNCIONES
inline double H(double a){							//	funcion de Hubble normalizada
return sqrt(Og/(a*a*a*a)+Om/(a*a*a)+Ok/(a*a)+Ol);}

inline double Hp(double a){						//	dericada logaritmica de la funcion de Hubble normalizada
return -0.5*(4*Og+3*Om*a+2*Ok*a*a)/(a*(Og+Om*a+Ok*a*a+Ol*a*a*a*a));}

inline double mE(double a){						//	rho_m/(a E)^2
return Om/(a*(Og+Om*a+Ok*a*a+Ol*a*a*a*a));}

inline double rP(double a){						//	(rho+3P)/(a^2 rho)
return (2*Og+Om*a-2.*Ol*a*a*a*a)/(a*a*(Og+Om*a+Ok*a*a+Ol*a*a*a*a));}

inline double g(double kx, double ky, double kz, double a){	//	funcion de Green para solucionar Poisson
  return -3*Om/(8*a*(sin(kx/2)*sin(kx/2)+sin(ky/2)*sin(ky/2)+sin(kz/2)*sin(kz/2)));
}

inline double kg_2(double kx, double ky, double kz){	//	funcion de Green para solucionar Poisson
  return 1.0/(sin(kx)*sin(kx)+sin(ky)*sin(ky)+sin(kz)*sin(kz));
}

inline double k_2(double kx, double ky, double kz){	//	funcion de Green para solucionar Poisson
  return 1.0/(kx*kx+ky*ky+kz*kz);
}

inline double window2(double x){	//	CIC fourier
  if(x!=0.)
    return pow(x/sin(x),4);
  else
    return 1.;
}



inline double window(double x,double y,double z){	//	CIC fourier delta
  double aux=1.;
  if(x!=0.)
    aux*=pow(x/sin(x),4);
  if(y!=0.)
    aux*=pow(y/sin(y),4);
  if(z!=0.)
    aux*=pow(z/sin(z),4);
  return aux;
}

inline double f(double a){						//	H0/a_dot
  return pow((Og/a+Om+Ok*a+Ol*a*a*a)/a,-0.5);}

inline double f1(double y2){						//	primera derivada de la funcion de crecimiento
return y2;}

inline double f2(double a,double y1,double y2){				//	segunda derivada de la funcion de crecimiento
return 3*Om/(2*a)*y1/(Ol*a*a*a*a+Ok*a*a+Om*a+Og)+(-3+0.5*(4*Og+3*Om*a+2*Ok*a*a)/(Ol*a*a*a*a+Ok*a*a+Om*a+Og))/a*y2;}

void RungeKutta(double a, double &y1,double &y2,double h){
  double k1y1, k2y1, k3y1, k4y1;
  double k1y2, k2y2, k3y2, k4y2;

  k1y1=h*f1(y2);		k1y2=h*f2(a,y1,y2);
  k2y1=h*f1(y2+0.5*k1y2);  	k2y2=h*f2(a+0.5*h,y1+0.5*k1y1,y2+0.5*k1y2);
  k3y1=h*f1(y2+0.5*k2y2);  	k3y2=h*f2(a+0.5*h,y1+0.5*k2y1,y2+0.5*k2y2);
  k4y1=h*f1(y2+k3y2);      	k4y2=h*f2(a+h,y1+k3y1,y2+k3y2);

  y1+=(k1y1+k2y1*2+k3y1*2+k4y1)/6;			//	Error en y del orden da^5
  y2+=(k1y2+k2y2*2+k3y2*2+k4y2)/6;			//	Error en y del orden da^5
}


