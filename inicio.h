// Usando la aproximacionde Zeldovich, genera una distribucion inicial de partículas a partir del espectro lineal generado por el CAMB en z=49, usando FFTW


void inicio(void){

cout<<endl<<"computing initial conditions via Zeldovich Approximation"<<endl;
//////////////////////////////// Leyendo los datos del espectro de potencias interpolado del CAMB

	FILE * D;
	D=fopen(file_CAMB,"r");

	if (D == NULL) {
		printf("Unable to open file_CAMB = %s\n",file_CAMB);
		exit(0);
	}
	else
		printf("reading file_CAMB = %s\n", file_CAMB);

	Int i=0,j;
	Doub aux1,aux2;
	while(fscanf(D,"%lf   %lf\n",&aux1, &aux2)!=EOF) i++;// cuenta el numero de lineas de entrada
	VecDoub K(i);
	VecDoub P(i);
	rewind(D);
	for(j=0;j<i;j++){
		fscanf(D,"%lf   %lf\n",&aux1, &aux2);
		K[j]=aux1;
		P[j]=aux2;}
	fclose(D);

	Spline_interp pk(K,P);				//	interpola           pk.interp(0.002)


////////////////////////////////////////////  dando valores a la funcion en el espacio de fourier

	int k,l,i2,j2,k2,i3,j3,k3,ijk,ijk2;		//	contadores
	double s,s1,s2,s3,s4,s5,s6;

	FILE * IN[1];					//	n_k archivos de salida (uno por cada escala deseada)
	char nomAr[22];
			
	sx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np3);
	sy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np3);
	sz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * np3);
	dftx = fftw_plan_dft_3d(np, np, np, sx, sx, FFTW_FORWARD, FFTW_ESTIMATE);
	dfty = fftw_plan_dft_3d(np, np, np, sy, sy, FFTW_FORWARD, FFTW_ESTIMATE);
	dftz = fftw_plan_dft_3d(np, np, np, sz, sz, FFTW_FORWARD, FFTW_ESTIMATE);


		cout<<endl<<"computing initial conditions"<<endl;
		Normaldev normal1 (0.,1.,152407);	//	Numeros aleatorios con distribucion gausiana     mu=0     sigma=1.	la tercera entrada es la semilla
		Normaldev normal2 (0.,1.,309);

		for(k=0;k<np;k++){				//		asignando los valores iniciales
	  	for(j=0;j<np;j++){
	  	for(i=0;i<np;i++){
			if((i==0) || (i==np/2))	{i2=i;}	else	{i2=np-i;}
			if((j==0) || (j==np/2))	{j2=j;}	else	{j2=np-j;}
			if((k==0) || (k==np/2))	{k2=k;}	else	{k2=np-k;}
			if(i<i2)	{i3=i;}		else	{i3=i2;}
			if(j<j2)	{j3=j;}		else	{j3=j2;}
			if(k<k2)	{k3=k;}		else	{k3=k2;}
	  	    ijk=i+np*(j+k*np); 			ijk2=i2+np*(j2+k2*np);		// contador
	  	    s=sqrt( pk.interp(dpinp*sqrt(i3*i3+j3*j3+k3*k3)) )    *k_2(i3*dpinp,j3*dpinp,k3*dpinp) *0.5;
	  	    if(i*i+j*j+k*k==0){
				sx[ijk][0]=0.;		sx[ijk][1]=0.;
				sy[ijk][0]=0.;		sy[ijk][1]=0.;
				sz[ijk][0]=0.;		sz[ijk][1]=0.;}
	  	    else if((i==0 || i==np/2) && (j==0 || j==np/2) && (k==0 || k==np/2)){	
			      s1=normal1.dev()*s;
			      sx[ijk][0]=i3*s1;		sx[ijk2][0]=i3*s1;	sx[ijk][1]=0.;		sx[ijk2][1]=0.;
			      sy[ijk][0]=j3*s1;		sy[ijk2][0]=j3*s1;	sy[ijk][1]=0.;		sy[ijk2][1]=0.;
			      sz[ijk][0]=k3*s1;		sz[ijk2][0]=k3*s1;	sz[ijk][1]=0.;		sz[ijk2][1]=0.;}
	  	    else{
	  	 	  	s1=normal1.dev()*s;
				s2=normal2.dev()*s;
	sx[ijk][0]=i3*s1;		sx[ijk2][0]=i3*s1;		sx[ijk][1]=i3*s2;		sx[ijk2][1]=-i3*s2;
	sy[ijk][0]=j3*s1;		sy[ijk2][0]=j3*s1;		sy[ijk][1]=j3*s2;		sy[ijk2][1]=-j3*s2;
	sz[ijk][0]=k3*s1;		sz[ijk2][0]=k3*s1;		sz[ijk][1]=k3*s2;		sz[ijk2][1]=-k3*s2;}
		}}}


//tansformada de fourier de las tres coordenadas
		fftw_execute(dftx);
		fftw_execute(dfty);
		fftw_execute(dftz);

//imprimiendo las condiciones iniciales

		sprintf(nomAr,"s_Lbox%d_Np%d.dat",Lbox,np);
		IN[0]=fopen(nomAr,"w+");

		for(k=0;k<np;k++){
		  for(j=0;j<np;j++){
		    for(i=0;i<np;i++){
		      ijk=i+j*np+k*np*np;
		      fprintf(IN[0],"%.7lf %.7lf %.7lf\n",sx[ijk][0],sy[ijk][0],sz[ijk][0]);}}}
		fclose(IN[0]);

/* // imprime la parte imaginaria
		FILE * R;
		R=fopen("r.dat","w+");
		for(k=0;k<np;k++){
		  for(j=0;j<np;j++){
		    for(i=0;i<np;i++){
		      ijk=i+j*np+k*np*np;
		      fprintf(R,"%.15lf  %.15lf  %.15lf\n",sx[ijk][1],sy[ijk][1],sz[ijk][1]);}}}
		fclose(R);
*/
		double auxx=0.,auxy=0.,auxz=0.;
		for(i=0;i<np3;i++){
			auxx+=sx[i][0]*sx[i][0];
			auxy+=sy[i][0]*sy[i][0];
			auxz+=sz[i][0]*sz[i][0];}
		cout<<auxx<<"	"<<auxy<<"	"<<auxz<<endl;

		auxx=0.;	auxy=0.;	auxz=0.;
		for(i=0;i<np3;i++){
			auxx+=sx[i][0];
			auxy+=sy[i][0];
			auxz+=sz[i][0];}
		cout<<auxx<<"	"<<auxy<<"	"<<auxz<<endl;

		auxx=0.;	auxy=0.;	auxz=0.;
		for(i=0;i<np3;i++){
			auxx+=sx[i][1]*sx[i][1];
			auxy+=sy[i][1]*sy[i][1];
			auxz+=sz[i][1]*sz[i][1];}
		cout<<auxx<<"	"<<auxy<<"	"<<auxz<<endl;


	fftw_destroy_plan(dftx);	fftw_destroy_plan(dfty);	fftw_destroy_plan(dftz);
	fftw_free(sx);			fftw_free(sy);			fftw_free(sz);

	cout<<endl;

}



