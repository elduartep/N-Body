
void anda(double &D_i,double &Dp_i,double &D_antes,double &Dp_antes){

	int i,j,k,aux2,mue,j2,k2;
	int p, ii, jj, kk, iii, jjj, kkk, ijk;
	//los dobles y triples son i+-1, j+-1, k+-1.    ijk es un contador que barre el numero de celdas

	double a,h;
	double aux,elk;		
	double gxp, gyp, gzp;
	double dx, dy, dz, tx, ty, tz; //variables auxliares para calcular el CIC
	double ff,fff;//es la funcion f evaluada
	double sx,sy,sz;//auxiliares
  	float f1,f2,f3,f4,f5,f6;

	char nomA[22];
	char nomAr[22];
	FILE * IN[2];	//archivo de entrada, y salida en z=0
	FILE * NOM[5];//los demas archivos de salida

	vx = (double*) fftw_malloc(sizeof(double) * np3 );
	vy = (double*) fftw_malloc(sizeof(double) * np3 );
	vz = (double*) fftw_malloc(sizeof(double) * np3 );
	x = (double*) fftw_malloc(sizeof(double) * np3 );
	y = (double*) fftw_malloc(sizeof(double) * np3 );
	z = (double*) fftw_malloc(sizeof(double) * np3 );

////////////////////			creando el plan
	rho = (double*) fftw_malloc(sizeof(double) * n3 );
	delta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n2* (n/2+1));
	dft  = fftw_plan_dft_r2c_3d(n, n, n, rho, delta, FFTW_MEASURE);
	dfti = fftw_plan_dft_c2r_3d(n, n, n, delta, rho, FFTW_MEASURE);






///////////////////
	if(GADGET_IC==1){
/////////////////////////////         leyendo el archivo con las posiciones y velocidades iniciales
	IN[0]=fopen(file_IC,"r");
	if (IN[0] == NULL) {
		printf("Unable to open %s\n", nomA);
		exit(0);
	}
	printf("reading %s\n", nomA);
	ijk=0;

	//double con=1.*(np/n+n*gama);
	double cons=D_i;
  double consta=pow(ai-aa*0.5,2)*Dp_antes/(pow(ai,2)*Dp_i);
  cout<<"velocity conversion given GADGET units = "<<consta<<endl;
  consta*=pow(ai,1.5)*n/Lbox*0.01;
	for(k=0;k<np;k++){
	  for(j=0;j<np;j++){
	    for(i=0;i<np;i++){
        // asigna las posiciones iniciales usando las condiciones iniciales GADGET 2LPTic
	    	fscanf(IN[0],"%f %f %f %f %f %f\n",&f1,&f2,&f3,&f4,&f5,&f6);//escanea las posiciones iniciales
  			f1*=1.*n/Lbox;			//	pasando a unidades de maquina
		        f2*=1.*n/Lbox;
		        f3*=1.*n/Lbox;
				//	asigna las posiciones iniciales usando la aproximacion de Zeldovich
		    	x[ijk]=fmod(1.*f1+0.5*0,n);
		    	y[ijk]=fmod(1.*f2+0.5*0,n);
		    	z[ijk]=fmod(1.*f3+0.5*0,n);
				//	asignando los momentos iniciales
			f4*=consta;
			f5*=consta;
			f6*=consta;
		    	vx[ijk]=(double) f4;
		    	vy[ijk]=(double) f5;
		    	vz[ijk]=(double) f6;
			ijk++;


}}}
  fclose(IN[0]);}



/////////////////////////////	Zeldovich Approximation
else{
        double cons=D_i *pow(double(np),-1.5)*sqrt(4*pi)   *128.0/Lbox     *pow(double(Lbox)/double(np),-2.5);
	aux=pow(ai-aa*0.5,2)*Dp_antes *pow(double(np),-1.5)*sqrt(4*pi)   *128.0/Lbox     *pow(double(Lbox)/double(np),-2.5);

	sprintf(nomA,"s_Lbox%d_Np%d.dat",Lbox,np);
	IN[0]=fopen(nomA,"r");
	if (IN[0] == NULL) {
		printf("Unable to open %s\n", nomA);
		exit(0);
	}
	printf("reading %s\n", nomA);
	sprintf(nomA,"ic_Lbox%d_Np%d.dat",Lbox,np);
	IN[1]=fopen(nomA,"w+");
	if (IN[1] == NULL) {
		printf("Unable to open %s\n", nomA);
		exit(0);
	}
	printf("writting %s\n", nomA);
	ijk=0;

	for(k=0;k<np;k++){
	  for(j=0;j<np;j++){
	    for(i=0;i<np;i++){

				//	asigna las posiciones iniciales usando la aproximacion de Zeldovich
	    	fscanf(IN[0],"%le %le %le\n",&sx,&sy,&sz);//escanea las posiciones iniciales
		    x[ijk]=fmod(n+1.0*i*n/np-cons*sx,n);
		    y[ijk]=fmod(n+1.0*j*n/np-cons*sy,n);
		    z[ijk]=fmod(n+1.0*k*n/np-cons*sz,n);
				//	asignando los momentos iniciales
		    vx[ijk]=-aux*sx;
		    vy[ijk]=-aux*sy;
		    vz[ijk]=-aux*sz;
	    	fprintf(IN[1],"%le %le %le\n",x[ijk]*Lbox/n,y[ijk]*Lbox/n,z[ijk]*Lbox/n);//escanea las posiciones iniciales
			ijk++;
}}}
  fclose(IN[0]);
  fclose(IN[1]);}


cout<<"evolving"<<endl;




	a=ai;
	aux2=0;
	while(a<=af){//ciclo principal

	cout<<"a = "<<a<<endl;




//////////////////////////////						calculo de la densidad////
	#pragma omp parallel for
	for(p=0;p<n3;p++){//loop sobre las celdas para borrar la densidad antigua
	  rho[p]=0.;}

////	#pragma omp parallel for private(i,j,k,dx,dy,dz,tx,ty,tz,ii,jj,kk)
	for(p=0;p<np3;p++){//loop sobre las particulas para calcular la densidad de las celdas
	  	i=floor(x[p]);	  	dx=x[p]-1.*i;	  	tx=1.-dx;
		j=floor(y[p]);		dy=y[p]-1.*j;		ty=1.-dy;
		k=floor(z[p]);		dz=z[p]-1.*k;		tz=1.-dz;
		//condiciones de frontera periodicas
		ii=(i+1)%n;				jj=(j+1)%n;			kk=(k+1)%n;
//		if((i<0) || (i>n) || (ii<0) || (ii>n) || (j<0) || (j>n) || (jj<0) || (jj>n) || (k<0) || (k>n) || (kk<0) || (kk>n)) cout<<endl<<"AQUI SI HAY UN ERROR, HAY UNA PARTICULA FUERA"<<endl;	
		//incrementando la densidade de las celdas adyacentes	(solo las componentes reales)	
	  rho[i+n*j	+n2*k]	+=tx*ty*tz;		rho[ii+n*j	+n2*k]	+=dx*ty*tz;
	  rho[i+n*jj	+n2*k]	+=tx*dy*tz;		rho[ii+n*jj	+n2*k]	+=dx*dy*tz;
	  rho[i+n*j	+n2*kk]	+=tx*ty*dz;		rho[ii+n*j	+n2*kk]	+=dx*ty*dz;
	  rho[i+n*jj	+n2*kk]	+=tx*dy*dz;		rho[ii+n*jj	+n2*kk]	+=dx*dy*dz;
}


	#pragma omp parallel for
	for(p=0;p<n3;p++){//loop sobre las celdas para borrar la densidad antigua
	  rho[p]*=pow((1.*n)/np,3);
	  rho[p]-=1.;}


		//calcula la transformada de Fourier de la sobredensidad
		fftw_execute(dft);


		//despues de tener tranformada de Fourier de la densidad podemos solucionar la ecuacion de 
		//Poisson en cada punto de red
//ijk=0;
#pragma omp parallel for private(j,i,ijk,aux)
		for(k=0;k<n;k++){//loop sobre los numeros de onda
		  for(j=0;j<n;j++){////				
		    for(i=0;i<(n/2+1);i++){

	  	    ijk=i+(n/2+1)*(j+k*n);
	//if((j==0) || (j==n/2))	{j2=j;}	else	{j2=n-j;}		if(j<j2)	{j2=j;}		else	;
	//if((k==0) || (k==n/2))	{k2=k;}	else	{k2=n-k;}		if(k<k2)	{k2=k;}		else	;
	//elk=dpin*sqrt(k2*k2+j2*j2+i*i);

					elk=dpin*sqrt(k*k+j*j+i*i);
					aux=g(i*dpin,j*dpin,k*dpin,a);//*window2(elk*0.5);
					delta[ijk][0]*=aux;			
					delta[ijk][1]*=aux;			
//          ijk++;
					}}}
		//arreglando a mano la componente 000
		delta[0][0]=0.;	delta[0][1]=0.;


		//transformada inversa del potencial
			fftw_execute(dfti);

//ahora vamos a calcular la aceleracion de cada particula interpolando la aceleracion en los puntos de red
		//el factor n³ en el deno. es la normaliz. de la transf. inversa
		ff=-0.5*f(a)*aa/n3;
		fff=f(a+0.5*aa)*aa/((a+0.5*aa)*(a+0.5*aa));
	int il,jl,kl;




double ax=0.,ay=0.,az=0.;
#pragma omp parallel for private(i,j,k,dx,dy,dz,tx,ty,tz,ii,jj,kk,iii,jjj,kkk,il,jl,kl,gxp,gyp,gzp,aux)
		for(p=0;p<np3;p++){//loop sobre las particulas para calcular la aceleracion de cada una
	i=floor(x[p]);  dx=x[p]-1.*i;		tx=1.-dx;		ii=(i+1)%n;		iii=(i+2)%n;	il=(n+i-1)%n;
	j=floor(y[p]);	dy=y[p]-1.*j;		ty=1.-dy;		jj=(j+1)%n;		jjj=(j+2)%n;	jl=(n+j-1)%n;
	k=floor(z[p]);	dz=z[p]-1.*k;		tz=1.-dz;		kk=(k+1)%n;		kkk=(k+2)%n;	kl=(n+k-1)%n;
			//calculo de la aceleracion para cada particula en cada direccion
			gxp	=		(rho[ii	+n*j	+n2*k]	-rho[il	+n*j	+n2*k ])*tx*ty*tz	+
						(rho[iii+n*j	+n2*k]	-rho[i	+n*j	+n2*k ])*dx*ty*tz	+
						(rho[ii	+n*jj	+n2*k]	-rho[il	+n*jj	+n2*k ])*tx*dy*tz	+
						(rho[iii+n*jj	+n2*k]	-rho[i	+n*jj	+n2*k ])*dx*dy*tz	+
						(rho[ii	+n*j	+n2*kk]	-rho[il	+n*j	+n2*kk])*tx*ty*dz	+
						(rho[iii+n*j	+n2*kk]	-rho[i	+n*j	+n2*kk])*dx*ty*dz	+
						(rho[ii	+n*jj	+n2*kk]	-rho[il	+n*jj	+n2*kk])*tx*dy*dz	+
						(rho[iii+n*jj	+n2*kk]	-rho[i	+n*jj	+n2*kk])*dx*dy*dz;

			gyp	=		(rho[i	+n*jj	+n2*k]	-rho[i	+n*jl	+n2*k])	*tx*ty*tz	+
						(rho[ii	+n*jj	+n2*k]	-rho[ii	+n*jl	+n2*k])	*dx*ty*tz	+
						(rho[i	+n*jjj	+n2*k]	-rho[i	+n*j	+n2*k])	*tx*dy*tz	+
						(rho[ii	+n*jjj	+n2*k]	-rho[ii	+n*j	+n2*k])	*dx*dy*tz	+
						(rho[i	+n*jj	+n2*kk]	-rho[i	+n*jl	+n2*kk])*tx*ty*dz	+
						(rho[ii	+n*jj	+n2*kk]	-rho[ii	+n*jl	+n2*kk])*dx*ty*dz	+
						(rho[i	+n*jjj	+n2*kk]	-rho[i	+n*j	+n2*kk])*tx*dy*dz	+
						(rho[ii	+n*jjj	+n2*kk]	-rho[ii	+n*j	+n2*kk])*dx*dy*dz;

			gzp	=		(rho[i	+n*j	+n2*kk]	-rho[i	+n*j	+n2*kl])*tx*ty*tz	+
						(rho[ii	+n*j	+n2*kk]	-rho[ii	+n*j	+n2*kl])*dx*ty*tz	+
						(rho[i	+n*jj	+n2*kk]	-rho[i	+n*jj	+n2*kl])*tx*dy*tz	+
						(rho[ii	+n*jj	+n2*kk]	-rho[ii	+n*jj	+n2*kl])*dx*dy*tz	+
						(rho[i	+n*j	+n2*kkk]-rho[i	+n*j	+n2*k ])*tx*ty*dz	+
						(rho[ii	+n*j	+n2*kkk]-rho[ii	+n*j	+n2*k ])*dx*ty*dz	+
						(rho[i	+n*jj	+n2*kkk]-rho[i	+n*jj	+n2*k ])*tx*dy*dz	+
						(rho[ii	+n*jj	+n2*kkk]-rho[ii	+n*jj	+n2*k ])*dx*dy*dz;

ax+=abs(gxp);	ay+=abs(gyp);	az+=abs(gzp);
			//actualizando las velocidades
			vx[p]+=ff*gxp;		vy[p]+=ff*gyp;		vz[p]+=ff*gzp;
//actualizando las posiciones

			aux=fmod(x[p]+vx[p]*fff+1.*n,1.*n);	x[p]=aux;
			aux=fmod(y[p]+vy[p]*fff+1.*n,1.*n);	y[p]=aux;
			aux=fmod(z[p]+vz[p]*fff+1.*n,1.*n);	z[p]=aux;
}



    a+=aa;		}//cierra el ciclo principal



	cout<<"	a= "<<a<<endl;

	sprintf(nomA,"xyz.dat");
	IN[1]=fopen(nomA,"w+");
	for(p=0;p<np3;p++){fprintf(IN[1],"%.7le %.7le %.7le\n",x[p]*Lbox/n,y[p]*Lbox/n,z[p]*Lbox/n);}
	fclose(IN[1]);

	FILE * VEL;
	sprintf(nomA,"vel%c_%d_5.dat",muestra[mue],n_in);
	VEL=fopen(nomA,"w+");
	for(p=0;p<np3;p++){fprintf(VEL,"%.7le %.7le %.7le\n",vx[p]*Lbox/n,vy[p]*Lbox/n,vz[p]*Lbox/n);}
	fclose(VEL);
	




	fftw_free(x);	fftw_free(y);	fftw_free(z);		//posicion de las particulas
	fftw_free(vx);	fftw_free(vy);	fftw_free(vz);		//velocidad de las particulas
	fftw_free(rho);	fftw_free(delta);
	fftw_destroy_plan(dft);					fftw_destroy_plan(dfti);



	cout<<"acaba la evolucion "<<endl;

cout<<endl<<"evolution finished"<<endl;

}

